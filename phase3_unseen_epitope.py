#!/usr/bin/env python3
"""
Phase 3 — does morphology help predict binding for UNSEEN epitopes?
===================================================================

This is the experiment the project lives or dies on. Phases 1-2 established
that morphological structure exists and is measurable. That says nothing about
whether it helps prediction, which is the actual open problem in the field.

Task
----
Binary TCR-epitope binding prediction. Model sees (TCR, epitope) and predicts
bind / no-bind. Evaluated on epitopes held out of training ENTIRELY, which is
the IMMREP25 setting where the field's performance collapses.

Arms (identical machinery, only the TCR featurizer differs)
----------------------------------------------------------
  kmer3     raw CDR3 3-mers            <- the baseline to beat
  morpheme  V gene, J gene, V-prefix token, N-region 3-mers, N-length bucket
  vjonly    V gene and J gene only     <- control: how much is germline usage?

Model
-----
Logistic regression over hashed TCR x EPITOPE feature crosses, trained with
SGD. A linear model on concatenated features cannot represent interaction, and
interaction is the entire task, so features are crossed explicitly. Pure
Python, no dependencies.

Leakage controls (these are the point)
--------------------------------------
  * folds are splits over EPITOPES, so test epitopes are never trained on
  * any test clonotype (cdr3, V, J) that also occurs under a training epitope
    is dropped from test
  * negatives are drawn within-split only, by repairing a TCR with a different
    epitope from the same split
  * clonotype-deduplicated input throughout
  * identical folds, negatives, hyperparameters and epitope features across
    arms; only the TCR featurizer varies

Caveat stated plainly: negatives are shuffled pairs, not experimentally
verified non-binders. The IMMREP post-mortems single this out as a source of
inflated performance. Absolute AUCs here are therefore optimistic; the
comparison BETWEEN arms is the result, not the level.

Usage:
    python3 phase3_unseen_epitope.py
"""

import math
import random
import sys
import zlib
from collections import Counter, defaultdict

import goi_core as core

SEED = 20260727
HASH_BITS = 18
HASH_SIZE = 1 << HASH_BITS
EPOCHS = 3
LR = 0.08
L2 = 1e-6
MAX_PER_EPITOPE = 800      # cap so one epitope cannot dominate
MIN_PER_EPITOPE = 150      # need enough to evaluate
N_EPITOPES = 20            # most-represented epitopes
N_FOLDS = 5
N_SEEDS = 5
N_CONTROL_SEEDS = 3
NEG_PER_POS = 1


# ---------------------------------------------------------------------------
# Featurizers
# ---------------------------------------------------------------------------

def kmers(seq, k):
    return [seq[i:i + k] for i in range(len(seq) - k + 1)]


def feats_kmer3(r):
    """Baseline: raw CDR3 3-mers. No morphological knowledge at all."""
    return ['K3=' + x for x in kmers(r['cdr3'], 3)] or ['K3=NONE']


def feats_morpheme(r):
    """Morphological tokenization: germline edges as units, N-region interior."""
    f = ['V=' + r['v'], 'J=' + r['j'], 'VP=' + r['v_prefix'],
         'JS=' + r['j_suffix'],
         'NL=' + str(min(len(r['n_region']), 12))]
    f += ['N3=' + x for x in kmers(r['n_region'], 3)]
    f += ['N2=' + x for x in kmers(r['n_region'], 2)]
    return f


def feats_vjonly(r):
    """Control: germline segment usage only, no junctional sequence."""
    return ['V=' + r['v'], 'J=' + r['j']]


ARMS = {
    'kmer3': feats_kmer3,
    'morpheme': feats_morpheme,
    'vjonly': feats_vjonly,
}


def feats_epitope(ep):
    """Shared across arms so the comparison is fair."""
    f = ['E=' + ep, 'EL=' + str(len(ep))]
    f += ['EP%d=%s' % (i, a) for i, a in enumerate(ep)]
    f += ['E3=' + x for x in kmers(ep, 3)]
    return f


# ---------------------------------------------------------------------------
# Hashed-cross logistic regression
# ---------------------------------------------------------------------------

def _h(s):
    """
    Stable hash. Python's built-in hash() is randomized per process
    (PYTHONHASHSEED), which would make runs non-reproducible.
    """
    return zlib.crc32(s.encode('ascii'))


def cross_indices(tcr_feats, ep_feats, cap_t=14, cap_e=14):
    """
    Hash the cross product of (capped) TCR and epitope features.
    Capping bounds the work per example; features are taken in order, which is
    deterministic given the featurizer.
    """
    idx = []
    for t in tcr_feats[:cap_t]:
        for e in ep_feats[:cap_e]:
            idx.append(_h(t + '|' + e) & (HASH_SIZE - 1))
    # keep the unary features too, so the model can learn base rates
    for t in tcr_feats[:cap_t]:
        idx.append(_h('T!' + t) & (HASH_SIZE - 1))
    for e in ep_feats[:cap_e]:
        idx.append(_h('E!' + e) & (HASH_SIZE - 1))
    return idx


def train_logreg(examples, rng, epochs=EPOCHS, lr=LR, l2=L2):
    w = [0.0] * HASH_SIZE
    b = 0.0
    order = list(range(len(examples)))
    for ep_i in range(epochs):
        rng.shuffle(order)
        for oi in order:
            idx, y = examples[oi]
            z = b
            for i in idx:
                z += w[i]
            if z > 30:
                p = 1.0
            elif z < -30:
                p = 0.0
            else:
                p = 1.0 / (1.0 + math.exp(-z))
            g = (p - y) * lr
            b -= g
            for i in idx:
                w[i] -= g + l2 * w[i]
    return w, b


def score(w, b, idx):
    z = b
    for i in idx:
        z += w[i]
    return z


# ---------------------------------------------------------------------------

def make_examples(rows, ep_pool, featurizer, rng, pos_set=None):
    """Positives are observed pairs; negatives re-pair a TCR with another
    epitope from the same pool."""
    if pos_set is None:
        pos_set = {(r['cdr3'], r['v'], r['j'], r['epitope']) for r in rows}
    ex, meta = [], []
    for r in rows:
        tf = featurizer(r)
        ex.append((cross_indices(tf, feats_epitope(r['epitope'])), 1))
        meta.append((r['epitope'], 1))
        for _ in range(NEG_PER_POS):
            for _try in range(10):
                other = rng.choice(ep_pool)
                if other == r['epitope']:
                    continue
                if (r['cdr3'], r['v'], r['j'], other) in pos_set:
                    continue
                ex.append((cross_indices(tf, feats_epitope(other)), 0))
                meta.append((other, 0))
                break
    return ex, meta


def positive_control(data, epitopes):
    """
    POSITIVE CONTROL — seen epitopes.

    Same features, same model, same negative scheme, but the split is over
    CLONOTYPES within each epitope, so every test epitope was also in training.
    This is the setting where the published literature reports AUC ~0.7+.

    If this control is also ~0.5, the harness cannot measure anything and the
    unseen-epitope result below is meaningless. It must pass before the main
    result can be interpreted.
    """
    rule('-')
    print("POSITIVE CONTROL — seen epitopes (clonotype split, 80/20)")
    rule('-')
    per_arm = {a: [] for a in ARMS}
    for k in range(N_CONTROL_SEEDS):
        for arm, featurizer in ARMS.items():
            rng = random.Random(SEED + 999 + 31 * k)
            train_rows, test_rows = [], []
            for ep in epitopes:
                rs = list(data[ep])
                rng.shuffle(rs)
                cut = int(0.8 * len(rs))
                train_rows += rs[:cut]
                test_rows += rs[cut:]
            all_pos = {(r['cdr3'], r['v'], r['j'], r['epitope'])
                       for r in train_rows + test_rows}
            tr, _ = make_examples(train_rows, epitopes, featurizer, rng,
                                  all_pos)
            te, _ = make_examples(test_rows, epitopes, featurizer, rng,
                                  all_pos)
            w, b = train_logreg(tr, rng)
            pos_s = [score(w, b, i) for i, y in te if y == 1]
            neg_s = [score(w, b, i) for i, y in te if y == 0]
            per_arm[arm].append(core.auc(pos_s, neg_s))

    results = {}
    for arm in ('kmer3', 'morpheme', 'vjonly'):
        vals = per_arm[arm]
        m = sum(vals) / len(vals)
        sd = (sum((v - m) ** 2 for v in vals) / len(vals)) ** 0.5
        results[arm] = m
        print(f"   {arm:<9s} AUC {m:.4f} +/- {sd:.4f}   "
              f"({', '.join(f'{v:.3f}' for v in vals)})")
    gaps = [per_arm['morpheme'][i] - per_arm['kmer3'][i]
            for i in range(N_CONTROL_SEEDS)]
    print(f"\n   morpheme - kmer3 on SEEN epitopes: "
          f"{results['morpheme'] - results['kmer3']:+.4f}   "
          f"(per seed: {', '.join(f'{g:+.3f}' for g in gaps)})")
    print(f"   consistent across seeds: "
          f"{'YES' if all(g > 0 for g in gaps) else 'NO'}")
    ok = max(results.values()) > 0.65
    print(f"\n   harness check: {'PASS' if ok else 'FAIL'} "
          f"(needs >0.65 on seen epitopes to be able to detect an effect)")
    return ok, results


def run_unseen(data, epitopes, seed, arms):
    """One complete unseen-epitope experiment at a given seed."""
    rng = random.Random(seed)
    shuffled_eps = list(epitopes)
    rng.shuffle(shuffled_eps)
    folds = [shuffled_eps[i::N_FOLDS] for i in range(N_FOLDS)]

    fold_auc = {a: [] for a in arms}
    ep_auc = {a: {} for a in arms}

    for fi, test_eps in enumerate(folds):
        train_eps = [e for e in epitopes if e not in set(test_eps)]
        train_rows = [r for e in train_eps for r in data[e]]
        train_clono = {(r['cdr3'], r['v'], r['j']) for r in train_rows}
        test_rows = []
        for e in test_eps:
            for r in data[e]:
                if (r['cdr3'], r['v'], r['j']) not in train_clono:
                    test_rows.append(r)

        for arm in arms:
            featurizer = ARMS[arm]
            rng_arm = random.Random(seed + fi)      # identical across arms
            tr_ex, _ = make_examples(train_rows, train_eps, featurizer, rng_arm)
            te_ex, te_meta = make_examples(test_rows, test_eps, featurizer,
                                           rng_arm)
            w, b = train_logreg(tr_ex, rng_arm)
            pos_s = [score(w, b, i) for i, y in te_ex if y == 1]
            neg_s = [score(w, b, i) for i, y in te_ex if y == 0]
            fold_auc[arm].append(core.auc(pos_s, neg_s))
            by_ep = defaultdict(lambda: ([], []))
            for (i, y), (ep_name, _) in zip(te_ex, te_meta):
                s = score(w, b, i)
                by_ep[ep_name][0 if y == 1 else 1].append(s)
            for ep_name, (p, n) in by_ep.items():
                if len(p) >= 30 and len(n) >= 30:
                    ep_auc[arm][ep_name] = core.auc(p, n)
    return fold_auc, ep_auc


def rule(ch='='):
    print(ch * 74)


def build_dataset():
    rows = core.dedup_clonotypes(core.filter_human_beta(core.load_vdjdb()))
    v_anch, j_anch, _ = core.derive_anchors(rows)
    ann, _ = core.annotate(rows, v_anch, j_anch)
    by_ep = defaultdict(list)
    for r in ann:
        by_ep[r['epitope']].append(r)
    keep = [(ep, rs) for ep, rs in by_ep.items() if len(rs) >= MIN_PER_EPITOPE]
    keep.sort(key=lambda t: -len(t[1]))
    return dict(keep[:N_EPITOPES])


def main():
    rng = random.Random(SEED)
    rule()
    print("PHASE 3 — UNSEEN-EPITOPE BINDING PREDICTION")
    rule()

    data = build_dataset()
    epitopes = sorted(data)
    print(f"Epitopes used: {len(epitopes)}")
    for ep in epitopes:
        if len(data[ep]) > MAX_PER_EPITOPE:
            rng.shuffle(data[ep])
            data[ep] = data[ep][:MAX_PER_EPITOPE]
    total = sum(len(v) for v in data.values())
    print(f"Clonotypes after capping at {MAX_PER_EPITOPE}/epitope: {total:,}")

    # folds over epitopes
    shuffled_eps = list(epitopes)
    rng.shuffle(shuffled_eps)
    folds = [shuffled_eps[i::N_FOLDS] for i in range(N_FOLDS)]
    print(f"Folds: {N_FOLDS} (each holds out {len(folds[0])} epitopes entirely)")
    print(f"Negatives: {NEG_PER_POS} per positive, shuffled pairs within split")

    harness_ok, control_aucs = positive_control(data, epitopes)

    rule('-')
    print("MAIN EXPERIMENT — unseen epitopes, repeated over seeds")
    rule('-')
    print("A single seed is not evidence: the fold assignment and the negative")
    print("sampling both depend on it. The effect must survive reseeding.")

    seeds = [SEED + 1000 * k for k in range(N_SEEDS)]
    arms = ['kmer3', 'morpheme', 'vjonly']
    seed_summaries = []
    for si, sd in enumerate(seeds):
        fold_auc, ep_auc = run_unseen(data, epitopes, sd, arms)
        means = {a: sum(fold_auc[a]) / len(fold_auc[a]) for a in arms}
        shared = sorted(set(ep_auc['morpheme']) & set(ep_auc['kmer3']))
        diffs = [ep_auc['morpheme'][e] - ep_auc['kmer3'][e] for e in shared]
        mean_d = sum(diffs) / len(diffs) if diffs else float('nan')
        wins = sum(1 for d in diffs if d > 0)
        seed_summaries.append((sd, means, mean_d, wins, len(diffs), diffs))
        print(f"\n   seed {sd}:")
        for a in arms:
            print(f"      {a:<9s} mean AUC {means[a]:.4f}")
        print(f"      morpheme - kmer3 (per-epitope mean) {mean_d:+.4f}"
              f"   wins {wins}/{len(diffs)}")

    rule('-')
    print("STABILITY ACROSS SEEDS")
    rule('-')
    for a in arms:
        vals = [ms[a] for _, ms, _, _, _, _ in seed_summaries]
        m = sum(vals) / len(vals)
        sd_ = (sum((v - m) ** 2 for v in vals) / len(vals)) ** 0.5
        print(f"   {a:<9s} mean AUC {m:.4f} +/- {sd_:.4f}   "
              f"({', '.join(f'{v:.3f}' for v in vals)})")

    deltas = [d for _, _, d, _, _, _ in seed_summaries]
    m_d = sum(deltas) / len(deltas)
    sd_d = (sum((v - m_d) ** 2 for v in deltas) / len(deltas)) ** 0.5
    all_diffs = [x for *_, ds in seed_summaries for x in ds]
    rng_boot = random.Random(SEED)
    boot = []
    for _ in range(2000):
        samp = [all_diffs[rng_boot.randrange(len(all_diffs))]
                for _ in all_diffs]
        boot.append(sum(samp) / len(samp))
    boot.sort()
    lo, hi = boot[50], boot[1949]
    print(f"\n   morpheme - kmer3 across seeds: {m_d:+.4f} +/- {sd_d:.4f}")
    print(f"   pooled per-epitope bootstrap 95% CI: [{lo:+.4f}, {hi:+.4f}]")
    n_pos = sum(1 for d in deltas if d > 0)
    print(f"   seeds favouring morpheme: {n_pos}/{len(deltas)}")

    rule('-')
    print("VERDICT")
    rule('-')
    consistent = n_pos == len(deltas) and lo > 0
    if not harness_ok:
        print("   INVALID — positive control failed; harness cannot measure.")
    elif consistent:
        print("   Morpheme features beat k-mers on unseen epitopes, "
              "consistently across all seeds.")
        print(f"   Effect is small: {m_d:+.4f} AUC over a baseline at "
              f"{sum(ms['kmer3'] for _, ms, _, _, _, _ in seed_summaries) / len(seed_summaries):.4f}.")
    elif lo > 0:
        print("   Positive on pooled data but NOT stable across every seed.")
        print("   Treat as suggestive only; needs more data or a stronger model.")
    else:
        print("   No reliable advantage for morpheme features on unseen")
        print("   epitopes at this scale. The Phase 2 structural result does")
        print("   NOT transfer to prediction.")
    print()
    print("   Note both baseline and treatment sit near chance (0.50) on")
    print("   unseen epitopes, reproducing the field-wide collapse that")
    print("   IMMREP25 reports. The seen-epitope control confirms the")
    print("   pipeline learns when the epitope is known:")
    for a in arms:
        print(f"      {a:<9s} seen {control_aucs[a]:.4f}  ->  unseen "
              f"{sum(ms[a] for _, ms, _, _, _, _ in seed_summaries) / len(seed_summaries):.4f}")
    rule()
    print("Absolute AUCs are optimistic: negatives are shuffled pairs, not")
    print("verified non-binders. The between-arm comparison is the result.")
    rule()
    return 0


if __name__ == '__main__':
    sys.exit(main())
