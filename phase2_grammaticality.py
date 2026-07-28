#!/usr/bin/env python3
"""
Phase 2 — the falsification test
================================

Necessary condition for the whole project: if the N-region has phonotactic
structure, a learned model must be able to tell real N-regions from decoys
that preserve everything except residue order.

Decoy design (this is the whole experiment):

  shuffled  — the real N-region with its residues permuted. Preserves length
              AND amino-acid composition exactly. Only ORDER is destroyed.
              Beating this means order carries information.
  resampled — residues drawn i.i.d. from the global N-region background,
              matched for length. Destroys composition too. Easy negative;
              included only as a sanity ceiling.

Models compared:

  length    — length distribution only. Control; should be ~0.5 against
              shuffles, since shuffling preserves length. If it is not, the
              evaluation is broken.
  flat      — bigram model over N-regions, pooled across all V/J. Tests
              order without any morphological conditioning.
  morpheme  — bigram model conditioned on the (V, J) segment pair, backing
              off to the pooled model. Tests whether the morphological
              context adds information beyond flat order.

Discipline:
  - clonotype-deduplicated input
  - TRAIN AND TEST SPLIT BY STUDY, so a model cannot win by memorizing one
    lab's batch; the test studies are unseen sources
  - AUC via Mann-Whitney U

Usage:
    python3 phase2_grammaticality.py
"""

import math
import random
import sys
from collections import Counter, defaultdict

import goi_core as core

RANDOM_SEED = 20260727
MIN_N_LEN = 3          # need at least 2 bigrams to say anything about order
SMOOTH = 0.5           # add-k
INTERP = 0.7           # weight on the conditioned model in the backoff mix

AA = 'ACDEFGHIKLMNPQRSTVWY'


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------

class BigramLM:
    """Add-k smoothed bigram model over amino acids with ^ / $ boundaries."""

    def __init__(self, smooth=SMOOTH):
        self.counts = defaultdict(Counter)
        self.totals = Counter()
        self.smooth = smooth
        self.vocab = len(AA) + 1          # + '$'

    def add(self, seq):
        toks = '^' + seq + '$'
        for a, b in zip(toks, toks[1:]):
            self.counts[a][b] += 1
            self.totals[a] += 1

    def logprob(self, seq):
        toks = '^' + seq + '$'
        lp = 0.0
        for a, b in zip(toks, toks[1:]):
            num = self.counts[a][b] + self.smooth
            den = self.totals[a] + self.smooth * self.vocab
            lp += math.log(num / den)
        return lp


class ConditionedLM:
    """Per-(V,J) bigram model interpolated with a pooled backoff."""

    def __init__(self, interp=INTERP):
        self.pooled = BigramLM()
        self.by_pair = defaultdict(BigramLM)
        self.interp = interp

    def add(self, seq, v, j):
        self.pooled.add(seq)
        self.by_pair[(v, j)].add(seq)

    def logprob(self, seq, v, j):
        lp_pool = self.pooled.logprob(seq)
        sub = self.by_pair.get((v, j))
        if sub is None or sub.totals['^'] < 20:
            return lp_pool
        lp_cond = sub.logprob(seq)
        # interpolate in log space via log-sum-exp of weighted probabilities
        a = math.log(self.interp) + lp_cond
        b = math.log(1 - self.interp) + lp_pool
        m = max(a, b)
        return m + math.log(math.exp(a - m) + math.exp(b - m))


class LengthModel:
    def __init__(self):
        self.dist = Counter()
        self.total = 0

    def add(self, seq):
        self.dist[len(seq)] += 1
        self.total += 1

    def logprob(self, seq):
        return math.log((self.dist[len(seq)] + 0.5) / (self.total + 0.5 * 40))


# ---------------------------------------------------------------------------
# Decoys
# ---------------------------------------------------------------------------

def shuffled_decoy(seq, rng):
    """Permute residues. Same length, same composition, different order."""
    chars = list(seq)
    for _ in range(10):
        rng.shuffle(chars)
        cand = ''.join(chars)
        if cand != seq:
            return cand
    return None          # e.g. homopolymer; no valid decoy exists


def resampled_decoy(seq, background, rng):
    return ''.join(rng.choices(background[0], weights=background[1],
                               k=len(seq)))


# ---------------------------------------------------------------------------

def bootstrap_auc_gap(pos_a, neg_a, pos_b, neg_b, rng, n_boot=200):
    """
    Paired bootstrap over test items for AUC(a) - AUC(b).
    Returns (2.5th, 97.5th) percentile of the gap.
    """
    n = len(pos_a)
    gaps = []
    idx_range = range(n)
    for _ in range(n_boot):
        idx = [rng.randrange(n) for _ in idx_range]
        pa = [pos_a[i] for i in idx]
        na = [neg_a[i] for i in idx]
        pb = [pos_b[i] for i in idx]
        nb = [neg_b[i] for i in idx]
        gaps.append(core.auc(pa, na) - core.auc(pb, nb))
    gaps.sort()
    lo = gaps[int(0.025 * n_boot)]
    hi = gaps[min(n_boot - 1, int(0.975 * n_boot))]
    return lo, hi


def rule(ch='='):
    print(ch * 74)


def main():
    rng = random.Random(RANDOM_SEED)

    rule()
    print("PHASE 2 — GRAMMATICALITY (FALSIFICATION TEST)")
    rule()

    rows = core.dedup_clonotypes(core.filter_human_beta(core.load_vdjdb()))
    v_anch, j_anch, _ = core.derive_anchors(rows)
    ann, _ = core.annotate(rows, v_anch, j_anch)
    ann = [r for r in ann if len(r['n_region']) >= MIN_N_LEN]
    print(f"Annotated clonotypes with N-region >= {MIN_N_LEN} aa: {len(ann):,}")

    # ------------------------------------------------- split by study, not row
    study_of = {}
    for r in ann:
        study_of[id(r)] = sorted(r['studies'])[0]
    studies = sorted({s for s in study_of.values()})
    rng.shuffle(studies)
    n_test = max(1, len(studies) // 5)
    test_studies = set(studies[:n_test])
    train = [r for r in ann if study_of[id(r)] not in test_studies]
    test = [r for r in ann if study_of[id(r)] in test_studies]
    print(f"Studies: {len(studies)}  ->  train {len(studies) - n_test}, "
          f"test {n_test} (held out entirely)")
    print(f"Clonotypes: train {len(train):,}  test {len(test):,}")

    # ----------------------------------------------------------------- train
    length_m = LengthModel()
    flat_m = BigramLM()
    cond_m = ConditionedLM()
    for r in train:
        n = r['n_region']
        length_m.add(n)
        flat_m.add(n)
        cond_m.add(n, r['v'], r['j'])

    # Capacity control: same model class, same number of conditioning classes,
    # but the (V, J) labels are randomly permuted across clonotypes. If the
    # conditioned model's advantage were merely extra parameters, this control
    # would show the same gain. It should not.
    real_labels = [(r['v'], r['j']) for r in (train + test)]
    shuffled_labels = list(real_labels)
    rng.shuffle(shuffled_labels)
    fake_label = {}
    for r, lab in zip(train + test, shuffled_labels):
        fake_label[id(r)] = lab
    perm_m = ConditionedLM()
    for r in train:
        fv, fj = fake_label[id(r)]
        perm_m.add(r['n_region'], fv, fj)

    bg_counts = Counter()
    for r in train:
        bg_counts.update(r['n_region'])
    background = (list(bg_counts.keys()), list(bg_counts.values()))

    # ------------------------------------------------------------- evaluate
    results = {}
    for decoy_name in ('shuffled', 'resampled'):
        pos = {'length': [], 'flat': [], 'morpheme': [], 'permuted': []}
        neg = {'length': [], 'flat': [], 'morpheme': [], 'permuted': []}
        n_pairs = 0
        for r in test:
            real = r['n_region']
            if decoy_name == 'shuffled':
                dec = shuffled_decoy(real, rng)
            else:
                dec = resampled_decoy(real, background, rng)
            if not dec:
                continue
            n_pairs += 1
            pos['length'].append(length_m.logprob(real))
            neg['length'].append(length_m.logprob(dec))
            pos['flat'].append(flat_m.logprob(real))
            neg['flat'].append(flat_m.logprob(dec))
            pos['morpheme'].append(cond_m.logprob(real, r['v'], r['j']))
            neg['morpheme'].append(cond_m.logprob(dec, r['v'], r['j']))
            fv, fj = fake_label[id(r)]
            pos['permuted'].append(perm_m.logprob(real, fv, fj))
            neg['permuted'].append(perm_m.logprob(dec, fv, fj))

        rule('-')
        print(f"DECOY: {decoy_name}   ({n_pairs:,} real/decoy pairs, "
              f"held-out studies)")
        rule('-')
        notes = {
            'length': 'control — must be ~0.50 vs shuffles',
            'flat': 'pooled bigram model, no morphology',
            'morpheme': 'conditioned on real (V,J)',
            'permuted': 'CONTROL — conditioned on shuffled (V,J) labels',
        }
        print(f"   {'model':<10s} {'AUC':>7s}   interpretation")
        aucs = {}
        for m in ('length', 'flat', 'morpheme', 'permuted'):
            aucs[m] = core.auc(pos[m], neg[m])
            print(f"   {m:<10s} {aucs[m]:>7.4f}   {notes[m]}")
        results[decoy_name] = (aucs, pos, neg)

        # paired bootstrap on the AUC gap that matters
        lo, hi = bootstrap_auc_gap(pos['morpheme'], neg['morpheme'],
                                   pos['flat'], neg['flat'], rng)
        gap = aucs['morpheme'] - aucs['flat']
        print(f"\n   morpheme - flat = {gap:+.4f}   "
              f"95% bootstrap CI [{lo:+.4f}, {hi:+.4f}]"
              f"   {'SIGNIFICANT' if lo > 0 else 'not significant'}")
        cap = aucs['permuted'] - aucs['flat']
        print(f"   permuted - flat = {cap:+.4f}   "
              f"(how much of the gain is model capacity, not morphology)")

    rule()
    print("Reading the result:")
    print("  flat AUC >> 0.5 vs shuffled  => residue ORDER carries signal")
    print("  morpheme AUC > flat AUC      => (V,J) morphological context adds")
    print("                                 information beyond flat order")
    print("  morpheme ~= flat             => decomposition is cosmetic; the")
    print("                                 core hypothesis is NOT supported")
    rule()
    return 0


if __name__ == '__main__':
    sys.exit(main())
