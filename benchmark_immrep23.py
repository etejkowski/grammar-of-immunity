#!/usr/bin/env python3
"""
Benchmark on the official IMMREP23 challenge data
=================================================

Addresses the main reviewer objection to the phase1-3 results: the baseline was
internal, and evaluation used our own splits rather than a community benchmark.

Data (from https://github.com/justin-barton/IMMREP23):
  VDJdb_paired_chain.csv  official sample training set, positives only
  solutions.csv           official test set WITH labels and public/private flag

Official metric: Macro AUC0.1 — partial ROC AUC to FPR 0.1 with McClish
standardisation, computed per peptide, then averaged. Implemented via
sklearn roc_auc_score(max_fpr=0.1), which is McClish-standardised.

Negative generation follows the challenge protocol: TCRs are swapped between
peptides whose Levenshtein distance exceeds 3, at 5 negatives per positive.

Arms (identical model and epitope features; only TCR featurisation differs):
  kmer3      CDR3b 3-mers                       <- baseline to beat
  morpheme   V, J, V-prefix, J-suffix, N-region k-mers, N-length
  vjonly     V and J gene identity only         <- germline-usage control
  cdr123     all CDR loops, both chains         <- what strong methods use

Published-style reference baselines, implemented here rather than cited:
  tcrbase    nearest-neighbour CDR3b similarity to that peptide's known
             binders (the TCRbase approach used as an IMMREP baseline)
  tcrdist    TCRdist-style weighted CDR mismatch distance, nearest neighbour

Crucially, the test set contains 20 peptides of which 13 occur in the official
training data and 7 do not, giving a seen/unseen split within one benchmark.

Usage:
    .venv/bin/python benchmark_immrep23.py
"""

import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from sklearn.feature_extraction import FeatureHasher
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

sys.path.insert(0, '.')
import goi_core as core

SEED = 20260727
N_NEG = 5
HASH_FEATURES = 2 ** 18
TRAIN_CSV = 'benchmarks/immrep23_VDJdb_paired_chain.csv'
TEST_CSV = 'benchmarks/immrep23_solutions.csv'

BLOSUM_FALLBACK = 4.0


def rule(ch='='):
    print(ch * 74)


# ---------------------------------------------------------------------------
# Distances
# ---------------------------------------------------------------------------

def levenshtein(a, b):
    if a == b:
        return 0
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1,
                           prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def norm_sim(a, b):
    """Length-normalised similarity in [0,1] from edit distance."""
    if not a or not b:
        return 0.0
    return 1.0 - levenshtein(a, b) / max(len(a), len(b))


# ---------------------------------------------------------------------------
# Featurisers
# ---------------------------------------------------------------------------

def kmers(s, k):
    s = s or ''
    return [s[i:i + k] for i in range(len(s) - k + 1)]


def make_featurisers(v_anch, j_anch):
    def f_kmer3(r):
        return ['K3=' + x for x in kmers(r['CDR3b'], 3)] or ['K3=NONE']

    def f_morpheme(r):
        v = core.base_gene(str(r['Vb']))
        j = core.base_gene(str(r['Jb']))
        d = core.decompose(str(r['CDR3b']).upper(), v, j, v_anch, j_anch)
        if not d:
            # fall back to raw k-mers when the parse fails, and flag it
            return ['MORPH=UNPARSED', 'V=' + v, 'J=' + j] + \
                   ['K3=' + x for x in kmers(r['CDR3b'], 3)]
        vp, n, js = d
        f = ['V=' + v, 'J=' + j, 'VP=' + vp, 'JS=' + js,
             'NL=' + str(min(len(n), 12))]
        f += ['N3=' + x for x in kmers(n, 3)]
        f += ['N2=' + x for x in kmers(n, 2)]
        return f

    def f_vjonly(r):
        return ['V=' + core.base_gene(str(r['Vb'])),
                'J=' + core.base_gene(str(r['Jb']))]

    def f_cdr123(r):
        f = []
        for col, tag in (('CDR1a', 'A1'), ('CDR2a', 'A2'), ('CDR3a', 'A3'),
                         ('CDR1b', 'B1'), ('CDR2b', 'B2'), ('CDR3b', 'B3')):
            val = r.get(col)
            if isinstance(val, str) and val:
                f += [f'{tag}3=' + x for x in kmers(val, 3)]
        return f or ['CDR=NONE']

    return {'kmer3': f_kmer3, 'morpheme': f_morpheme,
            'vjonly': f_vjonly, 'cdr123': f_cdr123}


def f_epitope(r):
    pep = str(r['Peptide'])
    hla = str(r.get('HLA', ''))
    f = ['E=' + pep, 'EL=' + str(len(pep)), 'H=' + hla]
    f += ['EP%d=%s' % (i, a) for i, a in enumerate(pep)]
    f += ['E3=' + x for x in kmers(pep, 3)]
    return f


def cross(tf, ef, cap_t=16, cap_e=16):
    out = []
    for t in tf[:cap_t]:
        for e in ef[:cap_e]:
            out.append(t + '|' + e)
    out += ['T!' + t for t in tf[:cap_t]]
    out += ['E!' + e for e in ef[:cap_e]]
    return out


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

def build_training_pairs(train, rng):
    """
    Official train file has positives only. Generate negatives per the
    challenge protocol: swap TCRs between peptides with Levenshtein > 3,
    5 negatives per positive.
    """
    peps = sorted(train['Peptide'].unique())
    far = {}
    for p in peps:
        far[p] = [q for q in peps if q != p and levenshtein(p, q) > 3]

    rows, labels = [], []
    by_pep = {p: g.index.tolist() for p, g in train.groupby('Peptide')}
    for idx, r in train.iterrows():
        rows.append(r)
        labels.append(1)
        pool = far[r['Peptide']]
        if not pool:
            continue
        for _ in range(N_NEG):
            donor_pep = pool[rng.integers(len(pool))]
            cand = by_pep[donor_pep]
            donor = train.loc[cand[rng.integers(len(cand))]]
            fake = donor.copy()
            fake['Peptide'] = r['Peptide']
            fake['HLA'] = r['HLA']
            rows.append(fake)
            labels.append(0)
    return pd.DataFrame(rows).reset_index(drop=True), np.array(labels)


def macro_auc01(df, scores, label_col='Label'):
    """Official metric: per-peptide AUC0.1 (McClish), then arithmetic mean."""
    per = {}
    for pep, g in df.assign(_s=scores).groupby('Peptide'):
        y = g[label_col].values
        if len(np.unique(y)) < 2:
            continue
        per[pep] = roc_auc_score(y, g['_s'].values, max_fpr=0.1)
    return (float(np.mean(list(per.values()))) if per else float('nan')), per


# ---------------------------------------------------------------------------

def main():
    rng = np.random.default_rng(SEED)
    rule()
    print("IMMREP23 OFFICIAL BENCHMARK")
    rule()

    train = pd.read_csv(TRAIN_CSV)
    test = pd.read_csv(TEST_CSV)
    train = train[train['Target'] == 1].reset_index(drop=True)
    print(f"official train positives : {len(train):,} "
          f"({train['Peptide'].nunique()} peptides)")
    print(f"official test rows       : {len(test):,} "
          f"({test['Peptide'].nunique()} peptides, "
          f"{test['Label'].mean():.1%} positive)")

    seen_peps = set(train['Peptide'])
    test_peps = sorted(test['Peptide'].unique())
    unseen = [p for p in test_peps if p not in seen_peps]
    print(f"test peptides in training: {len(test_peps) - len(unseen)}")
    print(f"test peptides UNSEEN     : {len(unseen)}  {unseen}")

    # sanity: all-zero submission must score 0.5 per the README
    zero_macro, _ = macro_auc01(test, np.zeros(len(test)))
    print(f"\nmetric sanity check, all-zero predictions: "
          f"Macro AUC0.1 = {zero_macro:.4f} (README says 0.5)")

    print("\ngenerating negatives per challenge protocol "
          f"({N_NEG} per positive, Levenshtein > 3)...")
    tr_pairs, tr_y = build_training_pairs(train, rng)
    print(f"training pairs: {len(tr_pairs):,} "
          f"({tr_y.mean():.1%} positive)")

    # morphological anchors from the training TCRs themselves
    anchor_rows = [{'cdr3': str(c).upper(),
                    'v': core.base_gene(str(v)), 'j': core.base_gene(str(j))}
                   for c, v, j in zip(train['CDR3b'], train['Vb'], train['Jb'])
                   if isinstance(c, str)]
    v_anch, j_anch, _ = core.derive_anchors(anchor_rows, min_support=10)
    print(f"germline anchors derived from official train: "
          f"{len(v_anch)} V, {len(j_anch)} J")
    featurisers = make_featurisers(v_anch, j_anch)

    parsed = sum(1 for r in anchor_rows
                 if core.decompose(r['cdr3'], r['v'], r['j'], v_anch, j_anch))
    print(f"decomposable training CDR3b: {parsed:,}/{len(anchor_rows):,} "
          f"({100 * parsed / len(anchor_rows):.1f}%)")

    results = {}
    rule('-')
    print("LEARNED MODELS — logistic regression on hashed TCRxepitope crosses")
    rule('-')
    hasher = FeatureHasher(n_features=HASH_FEATURES, input_type='string',
                           alternate_sign=False)
    for arm, fz in featurisers.items():
        Xtr = hasher.transform(
            cross(fz(r), f_epitope(r)) for _, r in tr_pairs.iterrows())
        Xte = hasher.transform(
            cross(fz(r), f_epitope(r)) for _, r in test.iterrows())
        clf = LogisticRegression(max_iter=2000, C=1.0, solver='liblinear')
        clf.fit(Xtr, tr_y)
        s = clf.decision_function(Xte)
        macro, per = macro_auc01(test, s)
        results[arm] = (macro, per, s)
        print(f"   {arm:<9s} Macro AUC0.1 = {macro:.4f}")

    rule('-')
    print("REFERENCE BASELINES — published-style, no training")
    rule('-')
    # TCRbase-style: similarity to the known binders of that peptide
    binders = defaultdict(list)
    for c, p in zip(train['CDR3b'], train['Peptide']):
        if isinstance(c, str):
            binders[p].append(c)
    tb = []
    for _, r in test.iterrows():
        ref = binders.get(r['Peptide'])
        q = str(r['CDR3b'])
        tb.append(max((norm_sim(q, x) for x in ref), default=0.0) if ref else 0.0)
    macro_tb, per_tb = macro_auc01(test, np.array(tb))
    results['tcrbase'] = (macro_tb, per_tb, np.array(tb))
    print(f"   tcrbase   Macro AUC0.1 = {macro_tb:.4f}   "
          f"(nearest known binder, CDR3b)")

    # TCRdist-style: CDR1+CDR2+CDR3 both chains, CDR3 weighted x3
    ref_loops = defaultdict(list)
    cols = ['CDR1a', 'CDR2a', 'CDR3a', 'CDR1b', 'CDR2b', 'CDR3b']
    for _, r in train.iterrows():
        ref_loops[r['Peptide']].append([str(r.get(c, '')) for c in cols])
    td = []
    for _, r in test.iterrows():
        ref = ref_loops.get(r['Peptide'])
        if not ref:
            td.append(0.0)
            continue
        q = [str(r.get(c, '')) for c in cols]
        best = 0.0
        for cand in ref:
            sim = 0.0
            wsum = 0.0
            for qi, ci, w in zip(q, cand, (1, 1, 3, 1, 1, 3)):
                sim += w * norm_sim(qi, ci)
                wsum += w
            best = max(best, sim / wsum)
        td.append(best)
    macro_td, per_td = macro_auc01(test, np.array(td))
    results['tcrdist'] = (macro_td, per_td, np.array(td))
    print(f"   tcrdist   Macro AUC0.1 = {macro_td:.4f}   "
          f"(weighted CDR1/2/3 both chains)")

    # ------------------------------------------------------- seen vs unseen
    rule('-')
    print("SEEN vs UNSEEN PEPTIDES (the result that matters)")
    rule('-')
    print(f"   {'model':<10s} {'all':>8s} {'seen':>8s} {'unseen':>8s}")
    for arm in ('kmer3', 'morpheme', 'vjonly', 'cdr123', 'tcrbase', 'tcrdist'):
        macro, per, _ = results[arm]
        s_vals = [v for p, v in per.items() if p in seen_peps]
        u_vals = [v for p, v in per.items() if p not in seen_peps]
        s_m = np.mean(s_vals) if s_vals else float('nan')
        u_m = np.mean(u_vals) if u_vals else float('nan')
        print(f"   {arm:<10s} {macro:>8.4f} {s_m:>8.4f} {u_m:>8.4f}")

    print(f"\n   n peptides: seen "
          f"{len([p for p in test_peps if p in seen_peps])}, "
          f"unseen {len(unseen)}")

    # paired per-peptide comparison, morpheme vs kmer3
    _, per_m, _ = results['morpheme']
    _, per_k, _ = results['kmer3']
    _, per_v, _ = results['vjonly']
    rule('-')
    print("PAIRED PER-PEPTIDE COMPARISONS")
    rule('-')

    def paired(a, b, subset, name):
        shared = sorted(p for p in (set(a) & set(b)) if subset(p))
        if len(shared) < 3:
            print(f"   {name:<28s} too few peptides ({len(shared)})")
            return
        d = np.array([a[p] - b[p] for p in shared])
        boot = np.array([np.mean(rng.choice(d, size=len(d), replace=True))
                         for _ in range(5000)])
        lo, hi = np.percentile(boot, [2.5, 97.5])
        verdict = ('better' if lo > 0 else
                   'no sig. difference' if hi > 0 else 'worse')
        print(f"   {name:<28s} n={len(d):<3d} wins {int((d > 0).sum()):>2d}"
              f"  delta {d.mean():+.4f}  CI [{lo:+.4f}, {hi:+.4f}]  {verdict}")

    paired(per_m, per_k, lambda p: True, 'morpheme - kmer3, all')
    paired(per_m, per_k, lambda p: p in seen_peps,
           'morpheme - kmer3, seen only')
    paired(per_m, per_k, lambda p: p not in seen_peps,
           'morpheme - kmer3, UNSEEN only')
    paired(per_m, per_v, lambda p: True, 'morpheme - vjonly, all')
    paired(per_m, per_v, lambda p: p in seen_peps,
           'morpheme - vjonly, seen only')

    rule('-')
    print("NOTE on the reference baselines")
    rule('-')
    print("   tcrbase and tcrdist score exactly 0.5000 on unseen peptides")
    print("   because they have no reference binders for a peptide absent from")
    print("   training and emit a constant score. Nearest-neighbour methods are")
    print("   structurally incapable of unseen-epitope prediction; that is a")
    print("   property of the method class, not a tuning failure.")
    rule()
    return 0


if __name__ == '__main__':
    sys.exit(main())
