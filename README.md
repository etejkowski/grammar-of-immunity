# The Grammar of Immunity

A linguistics-based approach to T-cell receptor binding prediction.

## What Is This?

This project applies **formal linguistics** (morphology, phonotactics, grammar formalisms) to the problem of predicting what pathogens a T-cell receptor (TCR) will recognize. The core insight: TCR CDR3 sequences are not random strings — they have morphological structure (V-prefix + N-region + J-suffix), and the N-region is where binding specificity is decided.

The framing follows Jerne's 1984 Nobel lecture, *The Generative Grammar of the Immune System*, which proposed the linguistics/immunology analogy but never formalized it, and Vu et al. (2024, *Nature Computational Science*), which specified what a formalization would require but did not build one.

## Status

Phase 1 (dataset) and Phase 2 (falsification test) are done. Phase 3 (does this
improve binding prediction on unseen epitopes?) is the open question and the
only one that determines whether the project matters.

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
- `phase2_grammaticality.py` — the falsification test above
- `grammar_of_immunity_demo.py` — original illustrative demo (superseded by the phase scripts; kept for provenance)

## Quick Start

```bash
git clone https://github.com/antigenomics/vdjdb-db.git
python3 phase1_dataset.py
python3 phase2_grammaticality.py
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

## What would falsify this project

Phase 3: if morpheme-aware tokenization does not beat k-mer tokenization at
predicting binding for epitopes held out of training entirely, then the
grammar is real but useless for the problem that matters, and that negative
result should be published.
