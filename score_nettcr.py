#!/usr/bin/env python3
"""
Score the retrained NetTCR-2.2 pan model on the official IMMREP23 test set
==========================================================================

Completes the published-method comparison. prepare_nettcr_data.py built the
training split; NetTCR-2.2's own train_nettcr_2_2_pan.py produced the checkpoint
in nettcr_run/model/checkpoint/. This script loads that checkpoint, predicts on
the official test set, and scores it with the same Macro AUC0.1 implementation
benchmark_immrep23.py uses, so the number drops straight into Table 5.

Encoding replicates train_nettcr_2_2_pan.py exactly: BLOSUM50 20aa, right-padded
to the published input widths, divided by 5.

Why retrained rather than released weights: see prepare_nettcr_data.py. The
released weights are trained on IEDB/VDJdb/10X data that IMMREP23 also draws
from, so they would be evaluated partly on their own training data.

Usage:
    .venv/bin/python score_nettcr.py
"""

import json
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

sys.path.insert(0, 'NetTCR-2.2/src')
sys.path.insert(0, '.')

SEED = 20260727
MODEL = 'nettcr_run/model/checkpoint/pan_immrep23.h5'
TEST = 'nettcr_run/test.csv'
OFFICIAL_TEST = 'benchmarks/immrep23_solutions.csv'
OFFICIAL_TRAIN = 'benchmarks/immrep23_VDJdb_paired_chain.csv'
PER_PEPTIDE_JSON = 'immrep23_per_peptide.json'
OUT = 'nettcr_run/nettcr_predictions.csv'

# published input widths, train_nettcr_2_2_pan.py lines 90-96
WIDTHS = {'peptide': 12, 'A1': 7, 'A2': 8, 'A3': 22,
          'B1': 6, 'B2': 7, 'B3': 23}


def rule(ch='='):
    print(ch * 74)


def macro_auc01(peptides, labels, scores):
    """Official metric: per-peptide AUC0.1 (McClish), then arithmetic mean."""
    df = pd.DataFrame({'peptide': peptides, 'y': labels, 's': scores})
    per = {}
    for pep, g in df.groupby('peptide'):
        if len(np.unique(g['y'].values)) < 2:
            continue
        per[pep] = roc_auc_score(g['y'].values, g['s'].values, max_fpr=0.1)
    return (float(np.mean(list(per.values()))) if per else float('nan')), per


def main():
    if not os.path.exists(MODEL):
        print(f"missing checkpoint: {MODEL}")
        print("run prepare_nettcr_data.py, then NetTCR-2.2's "
              "train_nettcr_2_2_pan.py")
        return 1

    import tensorflow as tf
    from tensorflow import keras
    import keras_utils

    rule()
    print("NetTCR-2.2 (pan, retrained on IMMREP23) — official test set")
    rule()

    test = pd.read_csv(TEST)
    official = pd.read_csv(OFFICIAL_TEST)
    train = pd.read_csv(OFFICIAL_TRAIN)
    train = train[train['Target'] == 1]
    seen_peps = set(train['Peptide'])

    print(f"test rows scored          : {len(test):,} of {len(official):,} "
          f"official")
    if len(test) != len(official):
        print("  NOTE: rows differ; comparison to Table 5 is on the retained "
              "subset only")
    print(f"test peptides             : {test['peptide'].nunique()} "
          f"({test['binder'].mean():.1%} positive)")

    # --- encode exactly as training did -----------------------------------
    enc = keras_utils.blosum50_20aa

    def encode(col, width):
        return np.float32(
            keras_utils.enc_list_bl_max_len(test[col], enc, width) / 5)

    x = {'pep': encode('peptide', WIDTHS['peptide']),
         'a1': encode('A1', WIDTHS['A1']),
         'a2': encode('A2', WIDTHS['A2']),
         'a3': encode('A3', WIDTHS['A3']),
         'b1': encode('B1', WIDTHS['B1']),
         'b2': encode('B2', WIDTHS['B2']),
         'b3': encode('B3', WIDTHS['B3'])}

    # the checkpoint was saved with a custom metric; supply it to load
    def my_numpy_function(y_true, y_pred):
        try:
            return roc_auc_score(y_true, y_pred, max_fpr=0.1)
        except ValueError:
            return np.array([float(0)])

    def auc_01(y_true, y_pred):
        return tf.numpy_function(my_numpy_function, [y_true, y_pred],
                                 tf.float64)

    model = keras.models.load_model(MODEL,
                                    custom_objects={'auc_01': auc_01})
    scores = model.predict(x, batch_size=64, verbose=0).ravel()

    y = test['binder'].values
    peps = test['peptide'].values
    macro, per = macro_auc01(peps, y, scores)

    s_vals = [v for p, v in per.items() if p in seen_peps]
    u_vals = [v for p, v in per.items() if p not in seen_peps]

    rule('-')
    print(f"   {'model':<12s} {'all':>8s} {'seen':>8s} {'unseen':>8s}")
    print(f"   {'nettcr2.2':<12s} {macro:>8.4f} "
          f"{np.mean(s_vals):>8.4f} {np.mean(u_vals):>8.4f}")
    print(f"\n   n peptides: seen {len(s_vals)}, unseen {len(u_vals)}")
    rule('-')
    print("PER-PEPTIDE Macro AUC0.1")
    rule('-')
    for pep in sorted(per, key=lambda p: -per[p]):
        tag = 'seen' if pep in seen_peps else 'UNSEEN'
        print(f"   {pep:<14s} {per[pep]:.4f}  {tag}")

    pd.DataFrame({'peptide': peps, 'binder': y,
                  'nettcr_score': scores}).to_csv(OUT, index=False)
    print(f"\nwrote {OUT}")

    # ---------------- paired comparison against the other arms -------------
    if os.path.exists(PER_PEPTIDE_JSON):
        ref = json.load(open(PER_PEPTIDE_JSON))['per_peptide']
        rng = np.random.default_rng(SEED)
        rule('-')
        print("PAIRED PER-PEPTIDE COMPARISONS, 5,000 bootstrap resamples")
        rule('-')

        def paired(other, subset, name):
            shared = sorted(p for p in (set(per) & set(other)) if subset(p))
            if len(shared) < 3:
                print(f"   {name:<30s} too few peptides ({len(shared)})")
                return
            d = np.array([per[p] - other[p] for p in shared])
            boot = np.array([np.mean(rng.choice(d, size=len(d), replace=True))
                             for _ in range(5000)])
            lo, hi = np.percentile(boot, [2.5, 97.5])
            verdict = ('better' if lo > 0 else
                       'no sig. difference' if hi > 0 else 'worse')
            print(f"   {name:<30s} n={len(d):<3d} wins "
                  f"{int((d > 0).sum()):>2d}  delta {d.mean():+.4f}  "
                  f"CI [{lo:+.4f}, {hi:+.4f}]  {verdict}")

        for arm in ('kmer3', 'morpheme', 'cdr123', 'tcrdist'):
            if arm not in ref:
                continue
            paired(ref[arm], lambda p: True, f'nettcr - {arm}, all')
            paired(ref[arm], lambda p: p in seen_peps,
                   f'nettcr - {arm}, seen')
            paired(ref[arm], lambda p: p not in seen_peps,
                   f'nettcr - {arm}, UNSEEN')
    else:
        print(f"({PER_PEPTIDE_JSON} not found; run benchmark_immrep23.py "
              f"for paired comparisons)")
    rule()
    return 0


if __name__ == '__main__':
    sys.exit(main())
