# The Grammar of Immunity

A linguistics-based approach to T-cell receptor binding prediction.

## What Is This?

This project applies **formal linguistics** (morphology, phonotactics, grammar formalisms) to the problem of predicting what pathogens a T-cell receptor (TCR) will recognize. The core insight: TCR CDR3 sequences are not random strings — they have morphological structure (V-prefix + N-region + J-suffix) and the N-region contains epitope-specific motifs that are invisible to flat sequence analysis.

## Key Finding

Morphological decomposition of 203K+ TCR sequences reveals **massive discriminative signals**:

| Comparison | Top Discriminative Bigram | Enrichment |
|---|---|---|
| Flu vs EBV | `IR` | **71.8x** |
| Flu vs CMV | `IR` | 33.3x |
| EBV vs Flu | `CF` | 26.3x |

## Files

- `grammar-of-immunity-research.md` — Full research plan, literature review, data analysis, 90-day roadmap
- `grammar_of_immunity_demo.py` — Working demo: loads VDJdb data, performs morphological decomposition

## Quick Start

```bash
git clone https://github.com/antigenomics/vdjdb-db.git
python3 grammar_of_immunity_demo.py
