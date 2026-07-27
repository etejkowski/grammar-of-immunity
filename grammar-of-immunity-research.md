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
decomposition reveals **highly discriminative motifs** (20-33x enrichment) in the N-region
that differentiate epitope specificities. This signal is invisible to naive tokenization.

---

## The Problem

### What is a T-Cell Receptor?

T-cells are immune cells that kill infected/cancerous cells. Each T-cell has a unique
receptor (TCR) on its surface — a protein that recognizes a specific fragment (peptide/epitope)
of a pathogen, presented on an MHC molecule.

### How TCRs Are Generated: V(D)J Recombination

Each TCR is assembled from gene segments:

