#!/usr/bin/env python3
"""
Boundary-trim control for the grammaticality result
===================================================

The objection this addresses
----------------------------
Empirical V anchors stop extending when the modal residue frequency falls below
0.80 (goi_core.derive_anchors). If a segment's true germline contribution runs
longer than the derived anchor, the leftover germline residues remain at the 5'
edge of what we call the N-region — and those residues are highly predictable
from V identity. A (V,J)-conditioned bigram model could therefore beat the flat
model partly by recognising boundary residues rather than junctional structure,
inflating the +0.1013 gap reported in Table 2.

The permuted-label control in phase2_grammaticality.py does NOT exclude this: it
destroys the V/J mapping, so it removes the leakage along with the signal.

The test
--------
Re-run the identical experiment with k residues removed from the N-region edges,
on the identical study split and seed, for four conditions:

    both k   trim k residues from each end
    v-side   trim k residues from the 5' (V) end only
    j-side   trim k residues from the 3' (J) end only

If the morpheme-minus-flat gap survives trimming, the effect is interior and the
Table 2 claim stands. If it collapses while the flat model's own order signal
holds up, the gap was boundary leakage.

Trimming also shortens sequences, which costs every model information, so the
flat AUC is reported alongside as the reference for how much loss is attributable
to length alone.

A direct leakage diagnostic is printed first: how predictable each N-region
position is from V (or J) identity, measured as the mean modal-residue frequency
within segment. Germline leakage shows up as high predictability at position 1
(or -1) relative to the interior.

Usage:
    python3 boundary_trim_control.py
"""

import random
import sys
from collections import Counter, defaultdict

import goi_core as core
from phase2_grammaticality import (RANDOM_SEED, MIN_N_LEN, BigramLM,
                                   ConditionedLM, shuffled_decoy,
                                   bootstrap_auc_gap)


def rule(ch='='):
    print(ch * 74)


def trim(seq, k_left, k_right):
    """Remove k_left residues from the 5' end and k_right from the 3' end."""
    end = len(seq) - k_right if k_right else len(seq)
    return seq[k_left:end]


def leakage_diagnostic(ann):
    """
    Mean modal-residue frequency at each N-region position, conditioned on the
    germline segment that abuts that end. High values near the edge mean the
    residue is largely determined by the segment, i.e. germline leakage.
    """
    rule('-')
    print("LEAKAGE DIAGNOSTIC — how segment-determined is each N-region "
          "position?")
    rule('-')
    print("   position   conditioned on   mean modal freq   n segments")
    for pos, side, key in ((0, 'V', 'v'), (1, 'V', 'v'), (2, 'V', 'v'),
                           (-1, 'J', 'j'), (-2, 'J', 'j'), (-3, 'J', 'j')):
        by_seg = defaultdict(Counter)
        for r in ann:
            n = r['n_region']
            if len(n) < 4:
                continue
            by_seg[r[key]][n[pos]] += 1
        fr, segs = [], 0
        for seg, c in by_seg.items():
            tot = sum(c.values())
            if tot < 50:
                continue
            fr.append(max(c.values()) / tot)
            segs += 1
        label = f"{pos:+d}" if pos < 0 else f"{pos + 1}"
        print(f"   {label:>8s}   {side:>14s}   {sum(fr) / len(fr):>15.3f}"
              f"   {segs:>10d}")
    print("\n   Interpretation: if position 1 is far more segment-determined")
    print("   than position 3, the 5' edge carries residual germline sequence.")


def run_condition(train, test, rng_seed, k_left, k_right, fake_label):
    """Train and evaluate flat / conditioned / permuted under one trim setting."""
    rng = random.Random(rng_seed)
    flat_m, cond_m, perm_m = BigramLM(), ConditionedLM(), ConditionedLM()

    n_train = 0
    for r in train:
        n = trim(r['n_region'], k_left, k_right)
        if len(n) < MIN_N_LEN:
            continue
        n_train += 1
        flat_m.add(n)
        cond_m.add(n, r['v'], r['j'])
        fv, fj = fake_label[id(r)]
        perm_m.add(n, fv, fj)

    pos = {'flat': [], 'morpheme': [], 'permuted': []}
    neg = {'flat': [], 'morpheme': [], 'permuted': []}
    for r in test:
        real = trim(r['n_region'], k_left, k_right)
        if len(real) < MIN_N_LEN:
            continue
        dec = shuffled_decoy(real, rng)
        if not dec:
            continue
        pos['flat'].append(flat_m.logprob(real))
        neg['flat'].append(flat_m.logprob(dec))
        pos['morpheme'].append(cond_m.logprob(real, r['v'], r['j']))
        neg['morpheme'].append(cond_m.logprob(dec, r['v'], r['j']))
        fv, fj = fake_label[id(r)]
        pos['permuted'].append(perm_m.logprob(real, fv, fj))
        neg['permuted'].append(perm_m.logprob(dec, fv, fj))

    aucs = {m: core.auc(pos[m], neg[m]) for m in pos}
    lo, hi = bootstrap_auc_gap(pos['morpheme'], neg['morpheme'],
                               pos['flat'], neg['flat'],
                               random.Random(rng_seed))
    return aucs, len(pos['flat']), n_train, (lo, hi)


def main():
    rng = random.Random(RANDOM_SEED)
    rule()
    print("BOUNDARY-TRIM CONTROL — is the (V,J) gain interior or edge leakage?")
    rule()

    rows = core.dedup_clonotypes(core.filter_human_beta(core.load_vdjdb()))
    v_anch, j_anch, _ = core.derive_anchors(rows)
    ann, _ = core.annotate(rows, v_anch, j_anch)
    ann = [r for r in ann if len(r['n_region']) >= MIN_N_LEN]
    print(f"Annotated clonotypes with N-region >= {MIN_N_LEN} aa: {len(ann):,}")

    leakage_diagnostic(ann)

    # identical split and label permutation to phase2_grammaticality.py
    study_of = {id(r): sorted(r['studies'])[0] for r in ann}
    studies = sorted({s for s in study_of.values()})
    rng.shuffle(studies)
    n_test = max(1, len(studies) // 5)
    test_studies = set(studies[:n_test])
    train = [r for r in ann if study_of[id(r)] not in test_studies]
    test = [r for r in ann if study_of[id(r)] in test_studies]

    real_labels = [(r['v'], r['j']) for r in (train + test)]
    shuffled_labels = list(real_labels)
    rng.shuffle(shuffled_labels)
    fake_label = {id(r): lab for r, lab in zip(train + test, shuffled_labels)}

    print(f"\nStudies: train {len(studies) - n_test}, test {n_test} "
          f"(held out entirely) — same split as Table 2")

    conditions = [('none (Table 2)', 0, 0),
                  ('both ends, 1', 1, 1),
                  ('both ends, 2', 2, 2),
                  ('both ends, 3', 3, 3),
                  ("5' (V) end, 1", 1, 0),
                  ("5' (V) end, 2", 2, 0),
                  ("3' (J) end, 1", 0, 1),
                  ("3' (J) end, 2", 0, 2)]

    rule('-')
    print("RESULTS — AUC vs order-shuffled decoys, held-out studies")
    rule('-')
    print(f"   {'trim':<16s} {'pairs':>7s} {'flat':>7s} {'morph':>7s} "
          f"{'perm':>7s} {'morph-flat':>11s} {'95% CI':>20s}")
    base_gap = None
    out = []
    for name, kl, kr in conditions:
        aucs, n_pairs, n_train, (lo, hi) = run_condition(
            train, test, RANDOM_SEED, kl, kr, fake_label)
        gap = aucs['morpheme'] - aucs['flat']
        if base_gap is None:
            base_gap = gap
        out.append((name, gap, aucs, n_pairs))
        print(f"   {name:<16s} {n_pairs:>7,d} {aucs['flat']:>7.4f} "
              f"{aucs['morpheme']:>7.4f} {aucs['permuted']:>7.4f} "
              f"{gap:>+11.4f}   [{lo:+.4f}, {hi:+.4f}]")

    rule('-')
    print("VERDICT")
    rule('-')
    print(f"   untrimmed morpheme-flat gap : {base_gap:+.4f}")
    for name, gap, aucs, _ in out[1:]:
        print(f"   {name:<16s} gap {gap:+.4f}  "
              f"({100 * gap / base_gap:.0f}% retained, "
              f"flat {aucs['flat']:.4f})")
    print("\n   If the gap retains most of its magnitude while flat AUC also")
    print("   falls, the loss is length/information, not edge leakage.")
    print("   If the gap collapses toward zero while flat AUC holds, the")
    print("   original gap was driven by residual germline residues.")
    rule()
    return 0


if __name__ == '__main__':
    sys.exit(main())
