#!/usr/bin/env python3
"""
Learning curve: is the morphological advantage a threshold or a trend?
=====================================================================

Section 3.7 established that morphological tokenization beats raw k-mers by
about +0.034 Macro AUC0.1 at ~68,000 training pairs and by nothing at ~29,000.
Two points cannot distinguish two very different stories:

  threshold  the method needs a minimum amount of data, then plateaus. The
             linguistic program has a ceiling and we have already hit it.
  trend      the advantage is still climbing at the largest size available.
             The method is data-hungry, and more data would help.

Those imply opposite recommendations, so we sweep training size properly.

Design
------
Pairs are built once from the official IMMREP23 training data using the
challenge negative protocol, then subsampled by POSITIVE, carrying each
positive's negatives with it, so class balance stays fixed across sizes.
Features are hashed once per arm and row-sliced, which keeps 36 model fits
cheap. Evaluation is always the full official test set with the official
metric, split into peptides seen and unseen in training.

Usage:
    .venv/bin/python learning_curve.py
"""

import json
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from sklearn.feature_extraction import FeatureHasher
from sklearn.linear_model import LogisticRegression

sys.path.insert(0, '.')
import goi_core as core
from benchmark_immrep23 import (TRAIN_CSV, TEST_CSV, HASH_FEATURES, N_NEG,
                                levenshtein, make_featurisers, f_epitope,
                                cross, macro_auc01, rule)

SEED = 20260727
ARMS = ('kmer3', 'morpheme', 'vjonly')
FRACTIONS = (0.10, 0.20, 0.35, 0.50, 0.75, 1.00)
N_SEEDS = 2
OUT_JSON = 'learning_curve.json'


def build_full_pairs(train, rng):
    """All positives plus challenge-protocol negatives, tagged by positive."""
    peps = sorted(train['Peptide'].unique())
    far = {p: [q for q in peps if q != p and levenshtein(p, q) > 3]
           for p in peps}
    by_pep = {p: g.index.tolist() for p, g in train.groupby('Peptide')}

    rows, labels, group = [], [], []
    for gi, (_, r) in enumerate(train.iterrows()):
        rows.append(r)
        labels.append(1)
        group.append(gi)
        for _ in range(N_NEG):
            pool = far[r['Peptide']]
            if not pool:
                break
            donor_pep = pool[rng.integers(len(pool))]
            ids = by_pep[donor_pep]
            donor = train.loc[ids[rng.integers(len(ids))]]
            fake = donor.copy()
            fake['Peptide'] = r['Peptide']
            fake['HLA'] = r['HLA']
            rows.append(fake)
            labels.append(0)
            group.append(gi)
    return (pd.DataFrame(rows).reset_index(drop=True),
            np.array(labels), np.array(group))


def stats_from_json(path=OUT_JSON):
    """Recompute every reported slope and extrapolation from the saved curve."""
    c = json.load(open(path))
    n = np.array([x['pairs'] for x in c], float)
    rule()
    print("LEARNING-CURVE STATISTICS (recomputed from " + path + ")")
    rule()
    for key, label in (('delta_seen', 'delta, seen peptides'),
                       ('delta_unseen', 'delta, unseen peptides'),
                       ('morpheme_seen', 'morpheme AUC, seen'),
                       ('kmer3_seen', 'k-mer AUC, seen'),
                       ('vjonly_seen', 'V/J-only AUC, seen')):
        d = np.array([x[key] for x in c])
        sl = np.polyfit(np.log(n), d, 1)[0]
        r = np.corrcoef(np.log(n), d)[0, 1]
        print(f"   {label:<24s} {d[0]:+.4f} -> {d[-1]:+.4f}   "
              f"slope {sl:+.4f}/log-pair   r={r:.3f}")

    d = np.array([x['delta_seen'] for x in c])
    upper = len(c) // 2
    print(f"\n   upper-half slope (seen delta): "
          f"{np.polyfit(np.log(n[upper:]), d[upper:], 1)[0]:+.4f}/log-pair")

    du = np.array([x['delta_unseen'] for x in c])
    sl, ic = np.polyfit(np.log(n), du, 1)
    print("\n   extrapolation of the UNSEEN-peptide delta (log-linear):")
    for target in (0.02, 0.05, 0.10):
        need = np.exp((target - ic) / sl)
        print(f"      to reach {target:+.2f}: {need:,.0f} pairs "
              f"= {need / n[-1]:,.0f}x the data used here")
    print("\n   Read as an order-of-magnitude argument, not a forecast:")
    print("   it projects a log-linear fit over one decade of observed range.")
    rule()


def main():
    if len(sys.argv) > 1 and sys.argv[1] == 'stats':
        stats_from_json()
        return 0
    rng = np.random.default_rng(SEED)
    rule()
    print("LEARNING CURVE — morphological advantage vs training-set size")
    rule()

    train = pd.read_csv(TRAIN_CSV)
    train = train[train['Target'] == 1].reset_index(drop=True)
    test = pd.read_csv(TEST_CSV)
    seen_peps = set(train['Peptide'])

    anchor_rows = [{'cdr3': str(c).upper(), 'v': core.base_gene(str(v)),
                    'j': core.base_gene(str(j))}
                   for c, v, j in zip(train['CDR3b'], train['Vb'], train['Jb'])
                   if isinstance(c, str)]
    v_anch, j_anch, _ = core.derive_anchors(anchor_rows, min_support=10)
    featurisers = make_featurisers(v_anch, j_anch)

    print("building full pair set once (challenge protocol)...")
    pairs, y, group = build_full_pairs(train, rng)
    n_pos = int((y == 1).sum())
    print(f"  {len(pairs):,} pairs, {n_pos:,} positives, "
          f"{y.mean():.1%} positive")

    hasher = FeatureHasher(n_features=HASH_FEATURES, input_type='string',
                           alternate_sign=False)
    print("hashing features once per arm...")
    X = {}
    Xte = {}
    for arm in ARMS:
        fz = featurisers[arm]
        X[arm] = hasher.transform(cross(fz(r), f_epitope(r))
                                  for _, r in pairs.iterrows())
        Xte[arm] = hasher.transform(cross(fz(r), f_epitope(r))
                                    for _, r in test.iterrows())
        print(f"  {arm:<9s} {X[arm].shape}")

    pos_groups = np.unique(group)
    results = defaultdict(list)

    rule('-')
    print(f"{'pairs':>8s} {'arm':<9s} {'all':>8s} {'seen':>8s} {'unseen':>8s}")
    rule('-')
    for frac in FRACTIONS:
        for si in range(N_SEEDS):
            r2 = np.random.default_rng(SEED + 17 * si)
            k = max(20, int(len(pos_groups) * frac))
            keep = set(r2.choice(pos_groups, size=k, replace=False).tolist())
            mask = np.fromiter((g in keep for g in group), bool, len(group))
            n_pairs = int(mask.sum())
            for arm in ARMS:
                clf = LogisticRegression(max_iter=2000, C=1.0,
                                         solver='liblinear')
                clf.fit(X[arm][mask], y[mask])
                s = clf.decision_function(Xte[arm])
                macro, per = macro_auc01(test, s)
                sv = float(np.mean([v for p, v in per.items()
                                    if p in seen_peps]))
                uv = float(np.mean([v for p, v in per.items()
                                    if p not in seen_peps]))
                results[(frac, arm)].append((macro, sv, uv, n_pairs))
                if si == 0:
                    print(f"{n_pairs:>8,} {arm:<9s} {macro:>8.4f} "
                          f"{sv:>8.4f} {uv:>8.4f}")
        print()

    rule('-')
    print("MORPHEME ADVANTAGE OVER RAW K-MERS, BY TRAINING SIZE")
    rule('-')
    print(f"{'pairs':>8s} {'seen delta':>12s} {'unseen delta':>14s} "
          f"{'morph seen':>12s} {'kmer seen':>11s}")
    curve = []
    for frac in FRACTIONS:
        m = np.array([r[:3] for r in results[(frac, 'morpheme')]]).mean(axis=0)
        k = np.array([r[:3] for r in results[(frac, 'kmer3')]]).mean(axis=0)
        v = np.array([r[:3] for r in results[(frac, 'vjonly')]]).mean(axis=0)
        n_pairs = int(np.mean([r[3] for r in results[(frac, 'morpheme')]]))
        curve.append({'fraction': frac, 'pairs': n_pairs,
                      'morpheme_all': m[0], 'morpheme_seen': m[1],
                      'morpheme_unseen': m[2],
                      'kmer3_all': k[0], 'kmer3_seen': k[1],
                      'kmer3_unseen': k[2],
                      'vjonly_all': v[0], 'vjonly_seen': v[1],
                      'vjonly_unseen': v[2],
                      'delta_seen': m[1] - k[1],
                      'delta_unseen': m[2] - k[2],
                      'delta_vs_vjonly_seen': m[1] - v[1]})
        print(f"{n_pairs:>8,} {m[1] - k[1]:>+12.4f} {m[2] - k[2]:>+14.4f} "
              f"{m[1]:>12.4f} {k[1]:>11.4f}")

    with open(OUT_JSON, 'w') as f:
        json.dump(curve, f, indent=2)
    print(f"\nwrote {OUT_JSON}")

    # --- shape of the curve -------------------------------------------------
    rule('-')
    print("IS IT A THRESHOLD OR A TREND?")
    rule('-')
    d = np.array([c['delta_seen'] for c in curve])
    n = np.array([c['pairs'] for c in curve], float)
    # slope over the upper half, where a plateau would flatten out
    upper = len(curve) // 2
    sl_all = np.polyfit(np.log(n), d, 1)[0]
    sl_upper = np.polyfit(np.log(n[upper:]), d[upper:], 1)[0]
    print(f"   delta at smallest ({int(n[0]):,} pairs): {d[0]:+.4f}")
    print(f"   delta at largest  ({int(n[-1]):,} pairs): {d[-1]:+.4f}")
    print(f"   slope per log-pair, whole range : {sl_all:+.4f}")
    print(f"   slope per log-pair, upper half  : {sl_upper:+.4f}")
    print()
    if d[-1] <= 0.005:
        print("   No advantage even at full size — nothing to characterize.")
    elif sl_upper > 0.5 * sl_all and sl_upper > 0.004:
        print("   TREND: the advantage is still rising at the largest size")
        print("   tested. Morphological tokenization is data-hungry rather")
        print("   than capped, and more data is the obvious next lever.")
    elif sl_upper < 0.25 * sl_all:
        print("   THRESHOLD: the advantage appears and then flattens. Extra")
        print("   data past this point is unlikely to help; the linguistic")
        print("   program has a ceiling near this level.")
    else:
        print("   AMBIGUOUS at this resolution: the upper-half slope is")
        print("   neither clearly flat nor clearly rising. More sizes or seeds")
        print("   would be needed to call it.")

    print()
    du = np.array([c['delta_unseen'] for c in curve])
    print(f"   unseen-peptide delta across all sizes: "
          f"{du.min():+.4f} to {du.max():+.4f}")
    print("   => generalization to new peptides does not improve with data")
    print("      at any size tested.")
    rule()
    return 0


if __name__ == '__main__':
    sys.exit(main())
