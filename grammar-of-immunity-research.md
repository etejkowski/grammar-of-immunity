# The Grammar of Immunity
## A Linguistics-Based Approach to TCR Binding Prediction

**Author**: Erick Tejkowski
**Date**: July 2026  
**Status**: Research Plan + Initial Data Validation

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [The Problem](#the-problem)
3. [Why Linguistics?](#why-linguistics)
4. [The Analogy (Formal)](#the-analogy-formal)
5. [Prior Work](#prior-work)
6. [Data Sources](#data-sources)
7. [Initial Data Analysis Results](#initial-data-analysis-results)
8. [The Research Plan (90 Days)](#the-research-plan-90-days)
9. [Tools & Infrastructure](#tools--infrastructure)
10. [Success Metrics](#success-metrics)
11. [Risks & Mitigations](#risks--mitigations)
12. [References](#references)

---

## Executive Summary

The adaptive immune system generates T-cell receptors (TCRs) through a stochastic
combinatorial process (V(D)J recombination) that produces astronomical diversity.
The central unsolved problem in computational immunology is: **given a TCR sequence,
predict what pathogen it will bind.**

Current approaches treat TCR sequences as generic protein strings and apply standard
NLP methods (transformers, BERT). They work for "seen" antigens but **fail completely
on novel pathogens** — exactly when prediction matters most (pandemics, cancer neoantigens).

**Our hypothesis**: TCR sequences have genuine linguistic structure (morphology, syntax,
semantics) that current models ignore. By applying formal linguistic analysis — treating
V(D)J recombination as a morphological derivation system and CDR3 motifs as semantic
units — we can build representations that generalize to unseen antigens.

**Initial validation**: Analysis of 203,308 TCR records from VDJdb shows that morphological
decomposition reveals **highly discriminative motifs** in the N-region that differentiate
epitope specificities. This signal is invisible to naive tokenization.

---

## The Problem

### What is a T-Cell Receptor?

T-cells are immune cells that kill infected/cancerous cells. Each T-cell has a unique
receptor (TCR) on its surface — a protein that recognizes a specific fragment (peptide/epitope)
of a pathogen, presented on an MHC molecule.

### How TCRs Are Generated: V(D)J Recombination

Each TCR is assembled from gene segments drawn from separate germline pools. For the
β chain: one V segment, one D segment, and one J segment are selected and joined. At each
junction the cell deletes a random number of nucleotides and inserts random ones (N/P
additions). The result is that the middle of the receptor — the part that actually touches
the antigen — is different in essentially every T-cell.

The **CDR3** (complementarity-determining region 3) spans this junction. It is the
principal determinant of binding specificity, and it is where the combinatorial diversity
concentrates: roughly 10¹⁵–10²⁰ theoretically possible sequences, of which any individual
carries perhaps 10⁷–10⁸.

### Why Prediction Is Hard

Two TCRs binding the same epitope may share little surface sequence similarity. Two nearly
identical TCRs may bind different epitopes. Known TCR–epitope pairs cover a tiny, heavily
biased slice of the space — a handful of epitopes (influenza M1, CMV pp65, EBV BMLF1)
account for a large fraction of all public data. So a model can score well by memorizing
the epitopes it has seen while learning nothing that transfers to a new one.

That failure mode is now well documented. Benchmarks find that models incorporating
multiple features outperform CDR3β-only models on seen epitopes, yet **all of them
struggle to generalize to unseen epitopes** [6].

---

## Why Linguistics?

The immune system is a generative system: a finite inventory of parts plus a set of
combination rules yielding unbounded, structured output. That is the same formal
situation natural language presents, and linguistics has spent a century building
machinery for exactly it.

Niels Jerne saw this and named it. His 1984 Nobel lecture was titled *The Generative
Grammar of the Immune System* [1], and it argued the analogy explicitly — but it stopped
at analogy. Forty years later the grammar has still not been written down.

The gap was restated precisely in 2024. Vu et al., a linguistics–immunology collaboration
at Oslo (the ImmunoLingo project), published a formalization of the "antibody language"
in *Nature Computational Science* [2]. Their argument: antibody language models borrow NLP
architectures but remain **domain-unspecific** — the tokenization is linguistically naive
and the training objective encodes no immune-specific grammar. They characterize what the
tokens and grammar would have to be, and they call for the implementation.

They did not build it. That is the opening this project takes.

The practical claim is narrower than the rhetoric: if CDR3s are morphologically derived
forms, then the correct unit of analysis is not the amino acid and not an arbitrary k-mer,
but the **morpheme** — the germline-contributed prefix, the stochastic junctional middle,
the germline-contributed suffix. Tokenizing at that boundary should expose regularities
that flat tokenization smears out.

---

## The Analogy (Formal)

| Linguistics | Immunology |
|---|---|
| Morphemes (word parts) | V, D, J gene segments |
| Word-formation rules (morphology) | V(D)J recombination |
| Words | CDR3 sequences |
| Sentences | Full receptor chains (α + β) |
| Meaning | Binding specificity |
| Grammar | **the unsolved part** |

### Level 1: Morphology (word formation)

V(D)J recombination behaves like derivational morphology:

- V segment = stem/root
- D segment = infix (β/heavy chains only)
- J segment = suffix
- N/P junctional additions = productive morphophonological processes

So a CDR3 is not "a sequence." It is a derived form with predictable germline-contributed
edges and a stochastic creative middle:

```
CDR3 = V-prefix + N-region + J-suffix
```

Current models tokenize by single amino acid or by fixed k-mer. Neither respects this
boundary.

### Level 2: Syntax (positional constraints)

CDR3s obey positional well-formedness conditions that behave syntactically:

- position 1 is essentially always C (cysteine) — like an obligatory capital letter
- the final position is F or W — like sentence-final punctuation
- certain residues are disfavored or forbidden at certain positions (= ungrammatical)
- permitted length ranges covary with V/J choice (= agreement)

### Level 3: Semantics (binding = meaning)

- CDR3s that bind the same epitope share motifs — synonyms sharing a root
- the same specificity is reachable by structurally dissimilar CDR3s — paraphrase
- some clonotypes are *public*, recurring across individuals — common vocabulary
- others are *private* to one person — idiolect

### Level 4: Pragmatics (context)

The MHC allele is the context of utterance. The same CDR3 can mean different things
against different MHC backgrounds — pragmatic ambiguity resolution.

---

## Prior Work

**Jerne (1984)** [1] — proposed the linguistics/immunology analogy in his Nobel lecture.
Programmatic, never formalized.

**Vu et al. (2024)** [2] — *Nature Computational Science*. Linguistics-based formalization
of the antibody language: defines tokens and grammar, diagnoses current antibody language
models as domain-unspecific, calls for grammar-aware implementations. A perspective paper;
no model.

**Lehner et al. (1995) / Moss, Callan et al. (1998)** [9] — the influenza M1
(GILGFVFTL) response was shown to be constrained to the BV17 family (TRBV19 in
current nomenclature) with a conserved **I/sRS(A)/S motif in CDR3β**, and the
structural basis of that conserved arginine-serine-serine was solved in 2003
[10]. This matters directly: any N-region analysis of GILGFVFTL that recovers
`IRSS` and its single-substitution variants is **reproducing a known result**.
Useful as pipeline validation; not a contribution.

**Xu et al. (2023), SPAN-TCR** [3] — *Cell Systems*. Scanning Parametrized by Normalized
TCR length. Entropic analysis identifies positional 2-mer motifs that **decrease the
entropy** of antigen-specific CDR3 groups, including motifs shared across different
specificities. This is the closest existing result to a proto-grammar, and it is
length-normalized rather than morphologically decomposed — it stops one step short.

**Dyrka, Pyzik, and colleagues** [4] — a sustained line of work fitting probabilistic
context-free grammars to protein sequences: binding-site detection, helix–helix contact
classification, amyloidogenic peptide recognition, and pCFG estimation under contact-map
constraints. Establishes that pCFGs are tractable and useful on protein sequence. Never
applied to the TCR specificity problem.

**IMMREP benchmarks** [5] — the community's yardstick. IMMREP25 focused exclusively on
unseen peptides: 126 named submissions predicting specificity of 1,000 TCRs against 20
unseen peptides restricted by HLA-A\*02:01 or HLA-B\*40:01. Independent assessments
concur that current methods do not generalize to truly unseen data [6].

**Where this project sits**: Jerne proposed the analogy, Vu et al. specified what a
formalization would require, SPAN-TCR found entropy-reducing positional motifs, and the
pCFG-for-proteins literature supplies the formal machinery. Nobody has connected the four.

---

## Data Sources

### VDJdb

- **What**: curated TCR sequences with known antigen specificity
- **Get it**: `git clone https://github.com/antigenomics/vdjdb-db`
- **Format**: TSV chunks, one file per publication, under `chunks/`
- **Fields used**: `cdr3.beta`, `v.beta`, `j.beta`, `antigen.epitope`, `mhc.a`, `species`, confidence score (0–3)

Field-to-linguistics mapping:

```
CDR3β:    CASSQDVGTGGVFALYF   <- the "word"
V gene:   TRBV14*01           <- prefix morpheme
J gene:   TRBJ1-5*01          <- suffix morpheme
Epitope:  GILGFVFTL           <- the "meaning"
MHC:      HLA-A*02:01         <- the context
```

### IEDB

- **What**: master database of known immune epitopes
- **API**: `https://query-api.iedb.org/tcell_search?linear_sequence=eq.SIINFEKL`
- **Tool**: TCRMatch, k-mer similarity matching over CDR3β

### IMGT

Germline reference for V/J segment sequences — required to do the morphological
decomposition correctly rather than by hand-coded table.

---

## Initial Data Analysis Results

Produced by `grammar_of_immunity_demo.py` in this repository. Stdlib Python, no
dependencies. **These figures are from a verified run, not from notes.**

### Corpus

| | |
|---|---|
| Total VDJdb records | 203,308 |
| Human TCRβ with full V/J/epitope annotation | 157,210 |

### Morphological decomposition

Decomposition uses a hand-built table of germline V-prefix and J-suffix contributions
(a stand-in for proper IMGT lookup, which is Phase 1 work).

| Epitope | Source | TCRs | Decomposed | Rate |
|---|---|---|---|---|
| GILGFVFTL | Influenza A (M1) | 14,152 | 13,357 | 94.4% |
| NLVPMVATV | CMV (pp65) | 18,925 | 16,709 | 88.3% |
| GLCTLVAML | EBV (BMLF1) | 8,954 | 6,016 | 67.2% |

The EBV rate is low because that set is dominated by V/J genes absent from the hand-built
table. Fixing this is exactly what IMGT integration buys.

### The finding: N-region motifs are epitope-specific

Influenza M1 (GILGFVFTL), top N-region motifs after stripping germline edges.
**These are database row counts, not clonotype counts** — `CASSIRSSYEQYF` alone
occupies 1,077 rows while being a single clonotype, so the percentages below
overstate motif prevalence severely. See the Phase 1 results below for
deduplicated figures.

```
'IRSS' — 1147 times (8.6%)
'IRST' —  317 times (2.4%)
'IRSA' —  292 times (2.2%)
'IRAS' —  164 times (1.2%)
'TRSS' —  139 times (1.0%)
```

Nearly one in ten flu-specific TCRs converges on the *same four-residue junctional
motif*, and the next four are single-substitution variants of it — `IRSS` → `IRST` →
`IRSA` → `IRAS`. That is a morphological paradigm: one root, systematic variation.

N-region length is also sharply constrained — 4 aa is the mode for flu (6,067 of 13,357),
while CMV peaks at 6 aa. Length covaries with specificity.

### Discriminative bigrams ("phonotactics")

Bigram enrichment, one epitope's N-regions against another's:

| Comparison | Top bigram | Enrichment |
|---|---|---|
| Flu vs EBV | `IR` | 62.5x |
| Flu vs CMV | `IR` | 24.5x |
| EBV vs Flu | `CF` | 30.8x |

### Caveats (read before quoting these numbers)

Two known issues, both to be addressed in Phase 1:

1. **Clonal redundancy.** The counts above include duplicate clonotype rows; VDJdb
   contains the same (CDR3, V, J) triple multiple times from different studies.
   Deduplicating drops `IR` in Flu vs EBV from **62.5x to 13.5x**. The signal survives,
   but raw-row figures are inflated roughly 5x. Deduplicated counts are the honest
   baseline and should be the ones reported in any writeup.

2. **Pseudocount inflation.** The enrichment ratio substitutes a 0.5 pseudocount for
   bigrams absent from the comparison set, which inflates rare bigrams. The EBV `CF`
   signal traces largely to ~95 sequences dominated by a single `ATDKLQLMKNCF` clonotype
   family. A minimum-count floor and a proper statistical test (Fisher exact with
   multiple-testing correction) are needed before any of these are called significant.

A third issue was found and fixed: the original decomposition mis-assigned the J-suffix
when a CDR3 carried an extra trailing residue, pushing the entire J region into the
N-region for 962 EBV sequences and manufacturing spurious `TN`/`NE`/`KN` "discriminative"
bigrams. Those artifacts are gone from the figures above.

**Interpretation.** The `IRSS` paradigm is real, survives deduplication, and
replicates across all 8 independent flu studies — but it was published in 1998
[9] and structurally explained in 2003 [10]. Recovering it validates the
pipeline; it is not a finding.

---

## Phase 1 and Phase 2 Results (completed)

`phase1_dataset.py` and `phase2_grammaticality.py` supersede the demo above.
Standard library only; deterministic under a fixed seed.

### Phase 1 — defensible dataset

| | |
|---|---|
| VDJdb rows | 203,308 |
| Human TCRβ, fully annotated | 157,206 |
| Unique clonotypes (CDR3, V, J, epitope) | 121,467 |
| Row inflation from duplication | 1.29x overall (2.07x for GILGFVFTL) |
| Decomposed | 115,951 (95.5%) |
| Dropped, malformed CDR3 (FR4 included) | 32 |
| J-suffix leakage into N-region | 0 (asserted) |

Germline anchors are now derived empirically — 71 V and 15 J segments, versus
the 40 hand-typed entries before and the 16 V / 13 J in VDJdb's bundled
reference — and validated against that reference where they overlap: 28
consistent, 1 discrepant (TRBV7-2, data `CASS` vs reference `CTSSL`).
Decomposition coverage for the EBV group rose from 67.2% to 91.3%.

### Phase 1 — batch confounding is the dominant methodological risk

| Epitope | Clonotypes | Studies | Largest single-study share |
|---|---|---|---|
| GILGFVFTL (flu) | 6,639 | 28 | 33.9% |
| NLVPMVATV (CMV) | 13,042 | 47 | 63.7% |
| GLCTLVAML (EBV) | 6,279 | 29 | **79.5%** (PMID 32184241) |

With Fisher exact + BH FDR + a count floor + cross-study replication, the flu
signals hold everywhere (`IR` 11.7x, q=5.8e-88, replicating in 8/8 studies with
≥200 clonotypes) while the EBV cysteine signals replicate in only 1 of 2
(`NC` 83.3x, 1/2; `CF` 18.0x, 1/2). That single dominant study has 10.8% of
N-regions containing cysteine against 0.4–0.6% in the flu and CMV groups. The
`CF` enrichment previously reported as a headline result is a batch artifact.

### Phase 2 — the falsification test passes

Can a model separate real N-regions from decoys that preserve length and
amino-acid composition and destroy **only order**? Trained on 320 studies,
evaluated on 79 entirely held-out studies, 16,304 real/decoy pairs:

| Model | AUC | Tests |
|---|---|---|
| length only | 0.5000 | control; must be 0.5 |
| flat bigram | 0.6149 | order, no morphology |
| **(V,J)-conditioned bigram** | **0.7162** | order + morphology |
| permuted-label control | 0.6115 | identical capacity, randomized labels |

morpheme − flat = **+0.1013**, 95% bootstrap CI [+0.0995, +0.1031].
permuted − flat = **−0.0035**, so none of the gain is attributable to the
conditioned model's extra parameters.

**What this establishes**: N-regions carry order-dependent structure, and
morphological context (the V/J pair) adds real information beyond flat
sequence statistics, generalizing to unseen studies. **What it does not
establish**: that any of this improves binding prediction. That is Phase 3.

### Phase 3 — the central hypothesis does not survive

Logistic regression over hashed TCR x epitope feature crosses (pure Python).
Identical folds, negatives, and hyperparameters across arms; only the TCR
featurizer differs. Negatives are shuffled pairs, so absolute levels are
optimistic and only between-arm comparisons are interpreted.

Seen epitopes (clonotype split, 3 seeds) — positive control:

| TCR features | AUC |
|---|---|
| raw CDR3 3-mers | 0.6589 +/- 0.0063 |
| morpheme | 0.7135 +/- 0.0017 |
| V/J genes only | 0.7089 +/- 0.0017 |

Unseen epitopes (5 seeds, epitopes withheld entirely, overlapping clonotypes
dropped):

| TCR features | AUC |
|---|---|
| raw CDR3 3-mers | 0.5052 +/- 0.0075 |
| morpheme | 0.5072 +/- 0.0065 |
| V/J genes only | 0.4991 +/- 0.0127 |

morpheme - kmer3 on unseen epitopes = **+0.0046 +/- 0.0104**, pooled
per-epitope bootstrap 95% CI **[-0.0084, +0.0175]**. Null. One seed produced
+0.0402 with a CI excluding zero; reseeding showed it was noise, which is why
the script sweeps seeds and reports stability.

Two conclusions, both negative and both useful:

1. **Morphology does not transfer.** The Phase 2 structural result is real but
   does not convert into unseen-epitope predictive power at this scale.
2. **The transferable part of the seen-epitope gain is germline usage, not
   grammar.** V/J identity alone reaches 0.7089 against morpheme's 0.7135, so
   the junctional N-region sequence — the "creative middle", the part the
   linguistic framing is really about — contributes roughly +0.005.

**Where the problem actually is**: the epitope side. No amount of TCR
tokenization tells a model how a *novel* peptide maps to the receptor features
it selects for. The next question worth asking is not "how do we tokenize
receptors" but "what representation of an epitope predicts which receptor
motifs it recruits."


---

## The Research Plan (90 Days)

### Phase 1: Foundation (Weeks 1–3)

Goal: data flowing, decomposition trustworthy.

- Clone and parse VDJdb into a working dataframe
- Exploratory analysis: CDR3 length by V/J, position-specific residue frequencies (= phonotactics), motifs within epitope groups
- **Replace the hand-built germline table with real IMGT lookup** so decomposition is principled and coverage is high across all V/J
- Deduplicate on (CDR3, V, J) and carry clone counts as weights rather than repeated rows
- Add a count floor and Fisher exact + FDR correction to the enrichment analysis

Deliverable: morphologically annotated, deduplicated dataset with defensible statistics.

### Phase 2: Grammar Induction (Weeks 4–7)

Goal: discover and formalize the rules.

- Per V–J combination, learn permitted length range, position-specific residue constraints, and bigram/trigram constraints
- **Falsification test**: can the learned grammar separate real CDR3s from shuffled/synthetic ones? If not, the grammar is vacuous and the project should stop here.
- Define a probabilistic context-free grammar:
  ```
  CDR3      -> V_prefix N_region J_suffix
  N_region  -> motif_1 | motif_2 | ...
  ```
  with probabilities estimated from data
- Compare against the existing protein pCFG literature [4]

Deliverable: a formal grammar generating valid CDR3s conditioned on V/J usage.

### Phase 3: Semantic Binding (Weeks 8–11)

Goal: does grammar-awareness actually predict better?

- **Baseline**: ESM-2 embeddings + classifier on raw CDR3
- **Model A**: same architecture, morphological tokenization (V-prefix | N-motif | J-suffix)
- **Model B**: Model A plus pCFG-derived features
- Evaluate on IMMREP-style splits: seen epitopes (should match state of the art) and **held-out unseen epitopes (the actual target)**
- Primary metric: AUC on unseen epitopes, with AUC₀.₁ for the high-precision regime

Deliverable: benchmarked model, honest comparison.

### Phase 4: Write Up (Week 12+)

- Works → submit to the next IMMREP round; write for *Bioinformatics* or *Nature Computational Science*
- Partially works → publish the grammar formalization as a methods contribution
- Fails → publish the negative result. "Linguistic structure does not transfer to
  specificity prediction" is a genuinely useful finding given how many people assume it does.

---

## Tools & Infrastructure

| Tool | Purpose |
|---|---|
| Python + pandas | data wrangling |
| IMGT / ANARCI | germline reference, IMGT numbering |
| OLGA | V(D)J recombination probability model |
| PyTorch or JAX | model training |
| ESM-2 | protein language model embeddings (baseline) |
| scirpy (Python) / immunarch (R) | repertoire analysis |
| NLTK | grammar/parsing utilities |

The current demo deliberately uses **only the standard library** so it runs anywhere with
Python 3.8+. That constraint ends at Phase 1.

---

## Success Metrics

Ordered by how much they'd matter:

1. **Primary**: AUC on unseen epitopes exceeds the published baselines. Anything above
   chance on genuinely held-out epitopes is a result, given that the field's honest
   answer today is "we can't."
2. **AUC₀.₁** — performance in the high-precision regime, which is what a clinical or
   vaccine-design use actually needs.
3. **Grammaticality discrimination**: the induced grammar separates real from shuffled
   CDR3s. A necessary condition; if this fails, nothing downstream is meaningful.
4. **Seen-epitope parity**: match, don't beat, existing methods. Losing here means the
   representation is throwing away information.
5. **Motif interpretability**: the grammar's high-probability productions correspond to
   motifs immunologists recognize. Buys credibility that a black box does not.

Explicit non-goal: leaderboard position. A small, well-characterized, honestly evaluated
result is worth more than an unreproducible one.

---

## Risks & Mitigations

| Risk | Mitigation |
|---|---|
| "I don't know enough immunology" | VDJdb docs plus lecture material gets you functional in about a week; the decomposition is already working |
| "The math is too hard" | pCFGs are tree-structured probability tables, not measure theory |
| "I can't compete with big labs" | They have compute; the linguistic framing is the asymmetry, and Vu et al. explicitly asked for it |
| Data sparsity for unseen epitopes | That *is* the problem statement — structure should generalize where memorization cannot |
| **Public data is clonally redundant and epitope-biased** | Deduplicate, weight by clone count, and never evaluate on an epitope that appears in training |
| **The grammar may be real but useless for prediction** | Phase 2's falsification test is the checkpoint; a negative result is publishable and honest |
| Hand-built germline tables silently corrupt decomposition | Already bit us once (J-suffix leakage). Now empirical anchors, validated against reference, with an assertion that no N-region contains a germline J suffix |
| **Epitope groups are dominated by single studies, so "specificity" signal is batch** | Report largest-study share per epitope; require cross-study replication before believing any enrichment; split train/test **by study** in every experiment |
| Rediscovering known motifs and mistaking them for findings | Literature-check every motif before claiming it (`IRSS` was published in 1998) |

---

## References

1. Jerne, N. K. (1985). The generative grammar of the immune system. Nobel lecture,
   8 December 1984. *Bioscience Reports* 5(6):439–451. Also *Science* (1985),
   doi:10.1126/science.4035345; PMID 3899210.
2. Vu, M. H., et al. (2024). Linguistics-based formalization of the antibody language as a
   basis for antibody language models. *Nature Computational Science*.
   doi:10.1038/s43588-024-00642-3. (ImmunoLingo project, University of Oslo.)
3. Xu, A. M., et al. (2023). Entropic analysis of antigen-specific CDR3 domains identifies
   essential binding motifs shared by CDR3s with different antigen specificities.
   *Cell Systems*. doi:10.1016/j.cels.2023.03.001; PMID 37001518. (SPAN-TCR)
4. Dyrka, W., Pyzik, M., et al. — probabilistic context-free grammars for protein
   sequence: binding-site detection (*BMC Systems Biology* 2007); helix–helix contact site
   classification and amyloidogenic peptide recognition (*Algorithms for Molecular Biology*
   2013, 8:31); pCFG estimation under contact-map constraints (*PeerJ* 2019, 7:e6559).
5. IMMREP25: Unseen Peptides — 126 named submissions, 1,000 TCRs, 20 unseen peptides
   restricted by HLA-A\*02:01 and HLA-B\*40:01.
6. Assessment of computational methods in predicting TCR–epitope binding recognition.
   *Nature Methods* (2025). doi:10.1038/s41592-025-02910-0. Finds that all evaluated
   methods struggle to generalize to unseen epitopes.
7. VDJdb — https://github.com/antigenomics/vdjdb-db
8. IEDB — https://www.iedb.org
9. A class I MHC-restricted recall response to a viral peptide is highly polyclonal
   despite stringent CDR3 selection. *J Immunol* (1998); PMID 9510187. Reports the
   BV17 (= TRBV19) constraint and the I/sRS(A)/S CDR3β motif for influenza M1
   GILGFVFTL — i.e. the motif this project's decomposition independently recovers.
10. A structural basis for immunodominant human T cell receptor recognition.
    *Nat Immunol* (2003); PMID 12796775. Vβ17Vα10.2 TCRs with a conserved
    arginine-serine-serine in CDR3β dominate the M1 response.
11. Tetramer-sort dataset contributing 79.5% of the GLCTLVAML clonotypes analyzed
    here — PMID 32184241 (VDJdb `reference.id` PMID:32184241, 23,666 rows,
    method.identification = tetramer-sort, no verification field).

---

## Provenance Note

Sections 3–12 of this document were reconstructed on 2026-07-27 after the original file was
truncated at 2,308 bytes during upload, losing everything after section 2. The
reconstruction draws on the originating session transcript for structure and argument.
All citations were independently re-verified against primary sources during
reconstruction. Section 7 was regenerated from a fresh verified run of
`grammar_of_immunity_demo.py` rather than transcribed, and therefore reports corrected
figures that differ from the pre-fix numbers.
