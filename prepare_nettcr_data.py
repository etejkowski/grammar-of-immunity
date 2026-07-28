#!/usr/bin/env python3
"""
Prepare IMMREP23 data in NetTCR-2.2 format
==========================================

Produces, in nettcr_run/:
    train.csv   official IMMREP23 training positives + challenge-protocol
                negatives, 80% of positives
    val.csv     the remaining 20%, for early stopping
    test.csv    the official IMMREP23 test set

Why retrain rather than use released weights
--------------------------------------------
NetTCR-2.2's published weights are trained on data compiled from IEDB, VDJdb and
10X Genomics, and its repository also distributes IMMREP 2022 benchmark training
data. IMMREP23's TCRs derive from those same sources, so scoring the IMMREP23
test set with released weights risks evaluating the model on its own training
data. Retraining on the identical split used by benchmark_immrep23.py is the only
comparison that means anything.

Column mapping: IMMREP23 CDR1a/CDR2a/CDR3a/CDR1b/CDR2b/CDR3b -> A1/A2/A3/B1/B2/B3.
Both use CDR3 positions 105-117, i.e. without the flanking C and F, so the
sequences transfer directly.

Negatives use the identical protocol and seed as benchmark_immrep23.py, so the
NetTCR arm sees exactly the training signal our own arms saw.

Usage:
    .venv/bin/python prepare_nettcr_data.py
"""

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, '.')
from benchmark_immrep23 import TRAIN_CSV, TEST_CSV, N_NEG, levenshtein

SEED = 20260727
OUTDIR = 'nettcr_run'
VAL_FRACTION = 0.20

COLMAP = {'CDR1a': 'A1', 'CDR2a': 'A2', 'CDR3a': 'A3',
          'CDR1b': 'B1', 'CDR2b': 'B2', 'CDR3b': 'B3'}

# NetTCR-2.2's published pan architecture fixes these input widths
# (src/train_nettcr_2_2_pan.py lines 90-96). Rows exceeding them are dropped
# rather than the architecture being widened, so the comparison is against the
# method as published.
CAPS = {'A1': 7, 'A2': 8, 'A3': 22, 'B1': 6, 'B2': 7, 'B3': 23,
        'peptide': 12}


def within_caps(df):
    ok = pd.Series(True, index=df.index)
    for col, cap in CAPS.items():
        ok &= df[col].str.len() <= cap
    return ok
NEEDED = list(COLMAP) + ['Peptide']


def to_nettcr(df, binder):
    out = df[NEEDED].rename(columns=COLMAP).copy()
    out['peptide'] = out.pop('Peptide')
    out['binder'] = binder
    return out


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    rng = np.random.default_rng(SEED)

    train = pd.read_csv(TRAIN_CSV)
    train = train[train['Target'] == 1].reset_index(drop=True)
    test = pd.read_csv(TEST_CSV)

    # drop rows lacking any CDR loop, which NetTCR cannot encode
    before = len(train)
    train = train.dropna(subset=NEEDED).reset_index(drop=True)
    print(f"training positives: {len(train):,} "
          f"(dropped {before - len(train)} with missing CDR loops)")

    peps = sorted(train['Peptide'].unique())
    far = {p: [q for q in peps if q != p and levenshtein(p, q) > 3]
           for p in peps}
    by_pep = {p: g.index.tolist() for p, g in train.groupby('Peptide')}

    rows, labels, group = [], [], []
    for gi, (_, r) in enumerate(train.iterrows()):
        rows.append(r)
        labels.append(1)
        group.append(gi)
        pool = far[r['Peptide']]
        if not pool:
            continue
        for _ in range(N_NEG):
            donor_pep = pool[rng.integers(len(pool))]
            ids = by_pep[donor_pep]
            donor = train.loc[ids[rng.integers(len(ids))]]
            fake = donor.copy()
            fake['Peptide'] = r['Peptide']
            rows.append(fake)
            labels.append(0)
            group.append(gi)

    pairs = pd.DataFrame(rows).reset_index(drop=True)
    y = np.array(labels)
    group = np.array(group)
    print(f"training pairs: {len(pairs):,} ({y.mean():.1%} positive)")

    # split by positive so a positive and its negatives stay together
    pos_ids = np.unique(group)
    rng.shuffle(pos_ids)
    n_val = int(len(pos_ids) * VAL_FRACTION)
    val_ids = set(pos_ids[:n_val].tolist())
    is_val = np.fromiter((g in val_ids for g in group), bool, len(group))

    def build(mask):
        sub = pairs[mask].reset_index(drop=True)
        lab = y[mask]
        out = sub[NEEDED].rename(columns=COLMAP)
        out = out.rename(columns={'Peptide': 'peptide'})
        out['binder'] = lab
        # NetTCR expects a partition column; a single partition is fine since
        # we supply explicit train and validation files.
        out['partition'] = 0
        out = out[['A1', 'A2', 'A3', 'B1', 'B2', 'B3', 'peptide',
                   'binder', 'partition']]
        keep = within_caps(out)
        if (~keep).any():
            print(f"    dropped {int((~keep).sum())} rows exceeding published "
                  f"input widths")
        return out[keep].reset_index(drop=True)

    tr = build(~is_val)
    va = build(is_val)
    te = test.dropna(subset=NEEDED).reset_index(drop=True)
    te_out = te[NEEDED].rename(columns=COLMAP).rename(
        columns={'Peptide': 'peptide'})
    te_out['binder'] = te['Label'].values
    te_out['partition'] = 0
    te_out = te_out[['A1', 'A2', 'A3', 'B1', 'B2', 'B3', 'peptide',
                     'binder', 'partition']]
    keep_te = within_caps(te_out)
    print(f"test rows exceeding published input widths: "
          f"{int((~keep_te).sum())}")
    te_out = te_out[keep_te].reset_index(drop=True)

    for name, df in (('train.csv', tr), ('val.csv', va),
                     ('test.csv', te_out)):
        path = os.path.join(OUTDIR, name)
        df.to_csv(path, index=False)
        print(f"  wrote {path}  {len(df):,} rows, "
              f"{df['binder'].mean():.1%} positive, "
              f"{df['peptide'].nunique()} peptides")

    # record which test rows were dropped, so scoring stays comparable
    dropped = len(test) - len(te_out)
    print(f"\ntest rows dropped for missing CDR loops: {dropped} "
          f"of {len(test):,}")
    print("Scoring will use only the retained rows for every arm compared.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
