# The Grammar of Immunity

A linguistics-based approach to T-cell receptor binding prediction.

## What Is This?

This project applies **formal linguistics** (morphology, phonotactics, grammar formalisms) to the problem of predicting what pathogens a T-cell receptor (TCR) will recognize. The core insight: TCR CDR3 sequences are not random strings — they have morphological structure (V-prefix + N-region + J-suffix) and the N-region contains epitope-specific motifs that are invisible to flat sequence analysis.

## Key Finding

Morphological decomposition of 203K+ TCR sequences reveals **massive discriminative signals**:

| Comparison | Top Discriminative Bigram | Enrichment |
|---|---|---|
| Flu vs EBV | `IR` | **62.5x** |
| Flu vs CMV | `IR` | 24.5x |
| EBV vs Flu | `CF` | 30.8x |

Figures are as reproduced by `grammar_of_immunity_demo.py` against VDJdb (203,308 records,
157,210 human TCRβ with full V/J/epitope annotation).

**Caveat, pending follow-up**: these counts include duplicate clonotype rows. Deduplicating
on `(cdr3, V, J)` reduces `IR` in Flu vs EBV from 62.5x to 13.5x. The signal survives, but the
raw-row figures are inflated roughly 5x by clonal redundancy, and the enrichment ratio uses a
0.5 pseudocount for absent bigrams, which amplifies rare ones. Both need addressing before
these numbers go into a writeup.

## Files

- `grammar-of-immunity-research.md` — Full research plan, literature review, data analysis, 90-day roadmap
- `grammar_of_immunity_demo.py` — Working demo: loads VDJdb data, performs morphological decomposition

## Quick Start

```bash
git clone https://github.com/antigenomics/vdjdb-db.git
python3 grammar_of_immunity_demo.py
