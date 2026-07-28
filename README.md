# The Grammar of Immunity

A linguistics-based approach to T-cell receptor binding prediction.

## What Is This?

This project applies **formal linguistics** (morphology, phonotactics, grammar formalisms) to the problem of predicting what pathogens a T-cell receptor (TCR) will recognize. The core insight: TCR CDR3 sequences are not random strings — they have morphological structure (V-prefix + N-region + J-suffix), and the N-region is where binding specificity is decided.

The framing follows Jerne's 1984 Nobel lecture, *The Generative Grammar of the Immune System*, which proposed the linguistics/immunology analogy but never formalized it, and Vu et al. (2024, *Nature Computational Science*), which specified what a formalization would require but did not build one.

## Status

Phases 1–3 are done. The headline: **the linguistic structure is real and
measurable, and it does not transfer to unseen-epitope prediction.**

## Result 4: Phase 3 — the hypothesis fails where it matters

`phase3_unseen_epitope.py`. Logistic regression over hashed TCR×epitope feature
crosses, pure Python. Identical folds, negatives, and hyperparameters across
arms; only the TCR featurizer differs.

**Seen epitopes** (clonotype split, 3 seeds) — positive control proving the
harness can measure:

| TCR features | AUC |
|---|---|
| raw CDR3 3-mers | 0.6589 ± 0.0063 |
| **morpheme (V, J, V-prefix, N-region k-mers)** | **0.7135 ± 0.0017** |
| V/J genes only | 0.7089 ± 0.0017 |

Morpheme beats raw k-mers by +0.0546, consistent across all 3 seeds. But note
the third row: **V and J gene identity alone gets 0.7089.** The entire
morphological advantage on seen epitopes is germline segment usage — a long
known predictor — and the junctional N-region sequence adds only +0.005 on top.

**Unseen epitopes** (5 seeds, epitopes held out entirely, overlapping
clonotypes dropped):

| TCR features | AUC |
|---|---|
| raw CDR3 3-mers | 0.5052 ± 0.0075 |
| morpheme | 0.5072 ± 0.0065 |
| V/J genes only | 0.4991 ± 0.0127 |

morpheme − kmer3 = **+0.0046 ± 0.0104**, pooled per-epitope bootstrap 95% CI
**[−0.0084, +0.0175]**. Null. A single seed gave +0.0402 with a CI excluding
zero; reseeding showed that was noise, which is why the script sweeps seeds.

Everything sits at chance, baseline included. This reproduces the field-wide
collapse IMMREP25 documents — the difference is that here it is reproduced with
a positive control demonstrating the same pipeline reaches 0.71 when the
epitope is known.

## Result 6: the morphological advantage is a data-scale effect

`negatives_robustness.py` retrains under four negative-generation schemes and
scores all of them on the same official IMMREP23 test set. Negative design is
the factor the IMMREP post-mortems blame for unstable results, so the
conclusions have to be shown invariant to it.

| Scheme | Decoy peptides | Pairs | morpheme − 3-mers (seen) | morpheme − V/J only |
|---|---|---|---|---|
| challenge | Levenshtein > 3 (official) | 67,872 | +0.0345 | +0.0013 |
| random | any other peptide | 67,872 | +0.0333 | −0.0023 |
| hard | Levenshtein ≤ 3 (adversarial) | 28,085 | −0.0022 | −0.0057 |
| matched | Levenshtein > 3, size-matched to `hard` | 29,957 | −0.0009 | −0.0068 |

The `hard` scheme serves negatives for only 3,729 of 11,312 positives, so it
changes training size and class balance as well as decoy difficulty. The
`matched` control holds size and balance fixed while restoring dissimilar
decoys, isolating the cause:

- `hard` (similar decoys, 28k pairs): **−0.0022**
- `matched` (dissimilar decoys, 30k pairs): **−0.0009**

Nearly identical, so **decoy similarity is not the driver**. What matters is
training set size: the morphological advantage over raw k-mers is +0.034 at
~68k training pairs and vanishes at ~29k, regardless of how hard the negatives
are. It is a data-scale effect, not a property of the representation, and
should not be reported as one.

Two conclusions survive every scheme:

1. **Nothing generalizes to unseen epitopes** — maximum across all arms and all
   four schemes is 0.5062.
2. **Morphological features are never meaningfully better than germline V/J
   identity alone** — +0.0013, −0.0023, −0.0057, −0.0068 across the four
   schemes.

## Result 5: confirmation on the official IMMREP23 benchmark

`benchmark_immrep23.py` repeats the test on community data — the official
IMMREP23 challenge training set and test set with released labels, scored with
the official metric (Macro AUC0.1, per-peptide partial ROC AUC to FPR 0.1 with
McClish standardisation). Negatives follow the challenge protocol (TCR swapping
between peptides at Levenshtein > 3, 5 per positive). The all-zero submission
scores exactly 0.5000, matching the challenge README, confirming the metric
implementation.

The test set contains 20 peptides, 13 present in the official training data and
7 absent — a seen/unseen split inside one benchmark.

| Model | All | Seen (13) | Unseen (7) |
|---|---|---|---|
| CDR3β 3-mers | 0.5602 | 0.6003 | 0.4858 |
| morpheme | 0.5893 | 0.6378 | 0.4992 |
| V/J genes only | 0.5802 | 0.6346 | 0.4791 |
| all CDR loops, both chains | 0.6039 | 0.6495 | 0.5192 |
| TCRbase-style nearest binder | 0.5919 | 0.6414 | 0.5000\* |
| TCRdist-style weighted CDR | **0.6281** | **0.6971** | 0.5000\* |

\* Nearest-neighbour methods emit a constant score for a peptide with no known
binders, so 0.5000 is structural: that method class cannot address unseen
epitopes at all.

Paired per-peptide comparisons:

| Comparison | n | Wins | Δ | 95% CI | |
|---|---|---|---|---|---|
| morpheme − 3-mers, all | 20 | 15 | +0.0291 | [+0.0093, +0.0505] | significant |
| morpheme − 3-mers, seen | 13 | 9 | +0.0375 | [+0.0086, +0.0677] | significant |
| morpheme − 3-mers, **unseen** | 7 | 6 | +0.0134 | [−0.0016, +0.0267] | **not significant** |
| morpheme − V/J only, all | 20 | 13 | +0.0091 | [−0.0166, +0.0338] | not significant |

Three independent confirmations of the internal results, now on community data
with a community metric:

1. Morphological tokenization beats raw k-mers (+0.0375 on seen peptides).
2. It does **not** beat V/J gene identity alone (+0.0091, CI crossing zero) —
   the germline-usage explanation holds.
3. On unseen peptides the advantage is not significant and every method sits
   near chance.

And one finding that goes against us: **our model is not state of the art.** A
TCRdist-style weighted CDR comparison using all six CDR loops of both chains
reaches 0.6281 overall and 0.6971 on seen peptides, well above the morpheme
model's 0.5893/0.6378. Using more of the receptor beats decomposing part of it
more cleverly.

### What this means

The bottleneck is not TCR tokenization. Better morphological representation of
the receptor cannot fix unseen-epitope generalization, because the failure is
on the **epitope side**: nothing in the model knows how a novel peptide maps to
the receptor features it should select for. Improving the TCR grammar is
optimizing the half of the problem that is not broken.

That is a useful negative result, and it redirects the research question from
"how do we tokenize TCRs?" to "what representation of an epitope predicts which
receptor motifs it will select?"

## Result 1: morphological structure is real and measurable

`phase2_grammaticality.py` asks the necessary question: can a model distinguish
real N-regions from decoys in which **only the residue order** has been
destroyed (same length, same amino-acid composition)?

Trained on 320 studies, evaluated on **79 entirely held-out studies**
(16,304 real/decoy pairs), clonotype-deduplicated:

| Model | AUC | What it tests |
|---|---|---|
| length only | 0.5000 | control — must be 0.5, and is exactly |
| flat bigram | 0.6149 | residue order, no morphology |
| **(V,J)-conditioned bigram** | **0.7162** | order + morphological context |
| permuted-label control | 0.6115 | same capacity, randomized (V,J) labels |

- Order carries information: 0.6149 vs 0.5 chance.
- Morphological conditioning adds **+0.1013 AUC** over flat, 95% bootstrap CI
  [+0.0995, +0.1031].
- The gain is **not** model capacity: shuffling the (V,J) labels while keeping
  the identical model class removes the entire advantage (−0.0035 vs flat).

So the morphological decomposition is not cosmetic. That is a necessary
condition for the project, not a sufficient one.

## Result 2: the pipeline reproduces known immunology (validation, not discovery)

Morphological decomposition of the influenza M1 epitope (GILGFVFTL) group
recovers a dominant N-region paradigm:

```
'IRSS' — 1147 clonotypes    'IRST' — 317    'IRSA' — 292    'IRAS' — 164
```

One root with systematic single-residue variation. **This is a rediscovery, not
a discovery.** The same motif was reported in 1998 — "the TCR beta-chain
repertoire ... constrained by the use of the BV17 family and the I/sRS(A)/S
amino acid motif in the CDR3 region" (PMID 9510187, where BV17 is TRBV19 in
current nomenclature) — and explained structurally in 2003 (PMID 12796775).

Recovering a 28-year-old published result from first principles is the best
available evidence that the decomposition captures real biology. It is not a
contribution.

## Result 3: most apparent "epitope-specific" signal is batch effect

With clonotype deduplication, Fisher exact tests, Benjamini-Hochberg FDR, a
count floor, and **cross-study replication**, the picture changes sharply:

| Comparison | Bigram | Ratio | q | Studies replicating |
|---|---|---|---|---|
| Flu > EBV | `IR` | 11.7x | 5.8e-88 | **8/8** |
| Flu > EBV | `MR` | 14.6x | 2.0e-21 | **8/8** |
| Flu > CMV | `RS` | 5.5x | 3.6e-178 | 8/8 |
| EBV > Flu | `NC` | 83.3x | 3.9e-25 | **1/2** |
| EBV > Flu | `CF` | 18.0x | 1.4e-05 | **1/2** |

The flu signals replicate across all eight independent studies with ≥200
clonotypes. The EBV cysteine signals do not: 79.5% of the GLCTLVAML group comes
from a single tetramer-sort study (PMID 32184241), which has 10.8% of N-regions
containing cysteine versus 0.4–0.6% elsewhere. Earlier versions of this README
reported `CF` at 26.3x as a headline finding. It is a batch artifact.

Run `python3 phase1_dataset.py` to regenerate; `enrichment.tsv` carries the
replication counts for every significant bigram.

## Files

- `grammar-of-immunity-research.md` — research plan, literature review, 90-day roadmap
- `goi_core.py` — data layer: loading, deduplication, empirical germline anchors, decomposition, Fisher/BH/AUC
- `phase1_dataset.py` — builds the annotated dataset, runs enrichment with FDR and replication testing
- `phase2_grammaticality.py` — the structural falsification test
- `phase3_unseen_epitope.py` — unseen-epitope prediction, with seen-epitope positive control and seed sweep
- `benchmark_immrep23.py` — official IMMREP23 benchmark, official metric, published-style reference baselines (requires numpy/scikit-learn/pandas)
- `negatives_robustness.py` — retrains under four negative-generation schemes, including a size-matched control
- `benchmarks/fetch.sh` — downloads the official IMMREP23 data
- `grammar_of_immunity_demo.py` — original illustrative demo (superseded by the phase scripts; kept for provenance)

## Quick Start

```bash
git clone https://github.com/antigenomics/vdjdb-db.git
python3 phase1_dataset.py
python3 phase2_grammaticality.py
python3 phase3_unseen_epitope.py    # ~5 min
```

For the official benchmark (the only part needing dependencies):

```bash
python3 -m venv .venv
.venv/bin/pip install numpy==1.26.4 scikit-learn==1.5.2 scipy==1.13.1 pandas==2.2.3
./benchmarks/fetch.sh
.venv/bin/python benchmark_immrep23.py    # ~7 min
```

Python 3.8+, standard library only. No dependencies.

## Method notes

Germline anchors are derived **empirically** rather than from a hand-typed
table: for each V segment, the 5' consensus is extended while the modal
residue holds ≥80% across clonotypes, and mirrored from the 3' end for J. This
covers 71 V and 15 J segments (vs 16 V / 13 J in VDJdb's bundled reference) and
is validated against that reference where it overlaps — 28 consistent, 1
discrepant (TRBV7-2: data says `CASS`, reference says `CTSSL`).

Records whose annotation extends past the CDR3 boundary (FR4 included, so the
germline J suffix appears twice) are detected and dropped rather than parsed —
32 such records exist.

## Where this stands

The project set out to test whether linguistic structure in TCR sequences would
generalize to unseen antigens. It ran the test and the answer, at this scale
with this model class, is no. The structure is real (Phase 2, +0.1013 AUC,
capacity-controlled) and it helps when the epitope is known (Phase 3 seen,
+0.0546 over k-mers), but it buys nothing for novel epitopes (+0.0046, CI
crossing zero).

Honest options from here, in descending order of how much they interest me:

1. **Move to the epitope side.** The failure is in mapping a novel peptide to
   the receptor features it selects. That is where the unsolved problem
   actually lives.
2. **Stronger model class.** A pCFG or a neural model might extract transferable
   structure that hashed feature crosses cannot. This is the Phase 2/3 plan in
   the research doc, and it is a real possibility — but note the baseline is at
   chance too, so the ceiling being tested is generalization itself, not
   representation.
3. **Publish the negative result.** "Morphological tokenization of TCR CDR3
   improves seen-epitope prediction but does not transfer to unseen epitopes,
   and the transferable component is germline V/J usage rather than junctional
   sequence" is a clean, useful claim with controls behind it.

Caveat on all absolute numbers: negatives are shuffled pairs, not verified
non-binders, which the IMMREP post-mortems flag as inflating performance. The
between-arm comparisons are the result; the levels are optimistic.
