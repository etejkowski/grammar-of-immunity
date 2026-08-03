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
collapse documented by the recent benchmarking assessments — Lu et al.
(*Nat Methods* 2026, 50 models across 21 datasets), Drost et al. (*Cell Genomics*
2025, 21 predictors via ePytope-TCR) and Liao et al. (arXiv:2606.04994), with the
IMMREP25 challenge reaching the same verdict — the difference is that here it is
reproduced with a positive control demonstrating the same pipeline reaches 0.71
when the epitope is known.

## Result 7: the size dependence is a trend, not a ceiling

`learning_curve.py` sweeps six training sizes (two seeds each) on the official
IMMREP23 data. The morphological advantage rises monotonically and never
flattens:

| Training pairs | morpheme (seen) | 3-mers (seen) | Δ seen | Δ unseen |
|---|---|---|---|---|
| 6,786 | 0.6108 | 0.6028 | +0.0080 | +0.0019 |
| 13,572 | 0.6211 | 0.6004 | +0.0207 | +0.0045 |
| 23,754 | 0.6290 | 0.6060 | +0.0230 | +0.0071 |
| 33,936 | 0.6414 | 0.6131 | +0.0283 | +0.0080 |
| 50,904 | 0.6365 | 0.6046 | +0.0318 | +0.0111 |
| 67,872 | 0.6378 | 0.6003 | **+0.0375** | +0.0134 |

+0.0118 Macro AUC0.1 per natural-log training pair (r = 0.985), with the
upper-half slope (+0.0131) slightly steeper than the whole-range slope. The
mechanism is visible when the arms are separated: **raw k-mers are flat across
the entire sweep** (0.6028 → 0.6003, slope +0.0009) while morphology climbs
(0.6108 → 0.6378, slope +0.0127). Morphological tokenization is not simply a
better representation — it is a more *data-responsive* one.

This reverses the pessimistic reading of Result 6. The ceiling we hit is a
property of how much annotated public data exists, not of the method.

One caveat keeps the result honest: the **V/J-only control rises fastest of all**
(+0.0167 per log-pair vs morphology's +0.0127). What extra data mainly buys is a
better estimate of germline segment preferences, not of junctional grammar — the
same conclusion Results 3 and 5 reached by different routes.

The unseen-peptide advantage also rises monotonically (+0.0048 per log pair,
r = 0.986), though it stays inside its confidence interval throughout. Taking the
slope at face value, a +0.05 unseen-peptide advantage would need roughly **2,400x
the available data** — which is the strongest argument in this repository for
moving to the antigen side rather than scaling the receptor side.

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
McClish standardization). Negatives follow the challenge protocol (TCR swapping
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
| NetTCR-2.2, retrained on this split | 0.5606 | 0.6003 | 0.4868 |

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
| NetTCR-2.2 − 3-mers, all | 20 | 12 | +0.0004 | [−0.0425, +0.0339] | not significant |
| NetTCR-2.2 − morpheme, all | 20 | 7 | −0.0287 | [−0.0783, +0.0104] | not significant |
| NetTCR-2.2 − TCRdist-style, all | 20 | 4 | −0.0675 | [−0.1164, −0.0272] | significant |

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

### A published deep model, retrained on the same split, does no better

The obvious objection to all of the above is model class: maybe hashed feature
crosses with logistic regression are simply too weak, and a real architecture
would generalize. So we retrained **NetTCR-2.2** (a CNN over all six CDR loops
of both chains, published defaults, 92,529 parameters) on the identical IMMREP23
training split and scored it with the identical metric.

```bash
.venv/bin/pip install tensorflow==2.15.1 matplotlib seaborn
git clone https://github.com/mnielLab/NetTCR-2.2.git   # third-party, own licence
.venv/bin/python prepare_nettcr_data.py
cd NetTCR-2.2/src && python train_nettcr_2_2_pan.py \
    --train_data ../../nettcr_run/train.csv \
    --val_data ../../nettcr_run/val.csv \
    --outdir ../../nettcr_run/model --model_name pan_immrep23
cd ../.. && .venv/bin/python score_nettcr.py
```

Result: **0.5606 overall, 0.6003 seen, 0.4868 unseen.** That is level with our
simplest arm, raw CDR3β 3-mers (+0.0004, CI [−0.0425, +0.0339]), below the
morpheme arm, and significantly below the TCRdist-style baseline. On all seven
unseen peptides it is at or below chance (best 0.5057).

The network is not undertrained: on its own held-out validation pairs, drawn from
the same peptides it trained on, it reaches AUC 0.9172 and AUC0.1 0.7821. It
fits the task well within distribution and collapses across peptides. The
unseen-epitope failure is therefore a property of the problem as currently posed,
not of our model class.

We retrained rather than using the released weights on purpose: those weights are
trained on IEDB/VDJdb/10X data that IMMREP23 also draws from, so scoring the
IMMREP23 test set with them risks evaluating the model on its own training data.

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

At the level of **distinct clonotypes**, influenza M1 (GILGFVFTL) specificity
recovers the published junctional signature:

| Epitope | N-region contains `RS` | Uses TRBV19 (=BV17) |
|---|---|---|
| Influenza M1 | **12.9%** | **41.0%** |
| CMV pp65 | 2.5% | 4.0% |
| EBV BMLF1 | 3.7% | 1.3% |

**This is a rediscovery, not a discovery.** The BV17 dominance of the HLA-A2
influenza M1 response was reported in 1991 and 1995 (PMID 1833769, PMID 7807026;
BV17 is TRBV19 in current nomenclature). A 1998 study then showed the response is
highly polyclonal *within* that constraint — 95 distinct CDR3β clonotypes sharing
the I/sRS(A)/S motif, with a power-law frequency distribution (PMID 9510187) —
and the structural basis was solved in 2003 (PMID 12796775). Recovering a result
of that age from first principles is evidence the decomposition tracks real
biology. It is not a contribution. Note that the 1998 polyclonality is exactly
what raw row counts destroy, which is the other half of the point here.

**Count rows and you will fool yourself.** The single sequence `CASSIRSSYEQYF`
appears in 1,077 VDJdb rows but is **one** clonotype. An earlier version of this
README reported row-based counts ("IRSS — 1147 times, 8.6%") as though they were
clonotype counts. After deduplication no single N-region exceeds nine
clonotypes; the motif is a distributed family with an `RS` core, not a handful of
frequent strings. See `fig6_motif_paradigm.png`.

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

## Figures

Generated by `make_figures.py` into `paper/figures/` (requires matplotlib):

| Figure | Shows |
|---|---|
| `fig1_decomposition.png` | what the morphological parse does, on real clonotypes |
| `fig2_grammaticality.png` | junctional structure is real, with capacity control |
| `fig3_seen_vs_unseen.png` | every method collapses on unseen peptides |
| `fig4_data_scale.png` | the advantage tracks data scale, not negative difficulty |
| `fig5_batch_effect.png` | one study can manufacture a "specificity" signal |
| `fig6_motif_paradigm.png` | rows vs clonotypes, and the 1998 motif recovered |
| `fig7_learning_curve.png` | the advantage grows with data and does not plateau |

## Files

- `paper/plain-language-guide.md` — **non-technical walkthrough of the study and all six figures** (start here)
- `paper/manuscript.md` — full manuscript draft
- `grammar-of-immunity-research.md` — research plan, literature review, 90-day roadmap
- `goi_core.py` — data layer: loading, deduplication, empirical germline anchors, decomposition, Fisher/BH/AUC
- `phase1_dataset.py` — builds the annotated dataset, runs enrichment with FDR and replication testing
- `phase2_grammaticality.py` — the structural falsification test
- `phase3_unseen_epitope.py` — unseen-epitope prediction, with seen-epitope positive control and seed sweep
- `benchmark_immrep23.py` — official IMMREP23 benchmark, official metric, published-style reference baselines (requires numpy/scikit-learn/pandas)
- `negatives_robustness.py` — retrains under four negative-generation schemes, including a size-matched control
- `learning_curve.py` — sweeps six training sizes to test whether the advantage is a threshold or a trend
- `prepare_nettcr_data.py` — builds the IMMREP23 split in NetTCR-2.2 format (identical protocol and seed)
- `score_nettcr.py` — scores the retrained NetTCR-2.2 checkpoint with the official metric, paired against every other arm
- `make_figures.py` — regenerates all seven figures (requires matplotlib)
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
   structure that hashed feature crosses cannot. Note that this is now partly
   answered: NetTCR-2.2, retrained on the identical split, lands at chance on
   unseen peptides too, so a CNN over all six loops is not the missing piece.
   The baseline is at chance as well, so the ceiling being tested is
   generalization itself, not representation.
3. **Publish the negative result.** "Morphological tokenization of TCR CDR3
   improves seen-epitope prediction but does not transfer to unseen epitopes,
   and the transferable component is germline V/J usage rather than junctional
   sequence" is a clean, useful claim with controls behind it.

Caveat on all absolute numbers: negatives are shuffled pairs, not verified
non-binders, which the IMMREP post-mortems flag as inflating performance. The
between-arm comparisons are the result; the levels are optimistic.
