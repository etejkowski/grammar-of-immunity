#!/usr/bin/env python3
"""
Grammar of Immunity — Morphological Decomposition Demo
======================================================

This script demonstrates the linguistic decomposition of T-cell receptor
CDR3 sequences into their morphological components:

    CDR3 = V-prefix + N-region + J-suffix

The N-region is the "semantic core" — the part that determines binding
specificity. By isolating it, we can discover epitope-specific motifs
that are invisible to flat sequence analysis.

Usage:
    git clone https://github.com/antigenomics/vdjdb-db.git
    python3 grammar_of_immunity_demo.py

Requirements:
    - Python 3.8+
    - No external packages needed (stdlib only)
"""

import os
import csv
from collections import Counter, defaultdict

# =============================================================================
# GERMLINE REFERENCE DATA
# =============================================================================

V_CONTRIBUTIONS = {
    'TRBV19': 'CASS', 'TRBV19*01': 'CASS', 'TRBV19-1': 'CASS',
    'TRBV19-1*01': 'CASS',
    'TRBV27': 'CASS', 'TRBV27*01': 'CASS',
    'TRBV20-1': 'CSAR', 'TRBV20-1*01': 'CSAR',
    'TRBV6-5': 'CASS', 'TRBV6-5*01': 'CASS',
    'TRBV6-1': 'CASS', 'TRBV6-1*01': 'CASS',
    'TRBV6-2': 'CASS', 'TRBV6-4': 'CASS', 'TRBV6-6': 'CASS',
    'TRBV28': 'CASS', 'TRBV28*01': 'CASS',
    'TRBV29-1': 'CSVG', 'TRBV29-1*01': 'CSVG',
    'TRBV9': 'CASS', 'TRBV9*01': 'CASS',
    'TRBV5-1': 'CASS', 'TRBV5-1*01': 'CASS',
    'TRBV5-6': 'CASS', 'TRBV5-8': 'CASS',
    'TRBV12-3': 'CASS', 'TRBV12-3*01': 'CASS',
    'TRBV12-4': 'CASS', 'TRBV12-4*01': 'CASS',
    'TRBV2': 'CASS', 'TRBV2*01': 'CASS',
    'TRBV7-9': 'CASS', 'TRBV7-9*01': 'CASS',
    'TRBV7-2': 'CASS', 'TRBV7-6': 'CASS', 'TRBV7-8': 'CASS',
    'TRBV14': 'CASS', 'TRBV14*01': 'CASS',
    'TRBV30': 'CAWS', 'TRBV30*01': 'CAWS',
    'TRBV4-1': 'CASS', 'TRBV4-2': 'CASS', 'TRBV4-3': 'CASS',
    'TRBV10-3': 'CASS', 'TRBV11-2': 'CASS',
    'TRBV13': 'CASS', 'TRBV15': 'CATS',
    'TRBV24-1': 'CATS', 'TRBV25-1': 'CASS',
}

J_CONTRIBUTIONS = {
    'TRBJ2-7': 'YEQYF', 'TRBJ2-7*01': 'YEQYF',
    'TRBJ2-1': 'NEQFF', 'TRBJ2-1*01': 'NEQFF',
    'TRBJ2-3': 'DTQYF', 'TRBJ2-3*01': 'DTQYF',
    'TRBJ1-2': 'YGYTF', 'TRBJ1-2*01': 'YGYTF',
    'TRBJ1-1': 'TEAFF', 'TRBJ1-1*01': 'TEAFF',
    'TRBJ2-2': 'GELFF', 'TRBJ2-2*01': 'GELFF',
    'TRBJ1-5': 'QPQHF', 'TRBJ1-5*01': 'QPQHF',
    'TRBJ2-5': 'ETQYF', 'TRBJ2-5*01': 'ETQYF',
    'TRBJ2-6': 'NVLTF', 'TRBJ2-6*01': 'NVLTF',
    'TRBJ1-3': 'NTIYF', 'TRBJ1-3*01': 'NTIYF',
    'TRBJ1-6': 'SPLYF', 'TRBJ1-6*01': 'SPLYF',
    'TRBJ1-4': 'NEKLF', 'TRBJ1-4*01': 'NEKLF',
    'TRBJ2-4': 'KNIQYF', 'TRBJ2-4*01': 'KNIQYF',
}

# =============================================================================
# MORPHOLOGICAL DECOMPOSITION
# =============================================================================

def decompose_cdr3(cdr3, v_gene, j_gene):
    """
    Morphological decomposition of a CDR3 sequence into:
        (V-prefix, N-region, J-suffix)
    """
    if not cdr3 or len(cdr3) < 5:
        return None
    v_prefix = None
    for vg in [v_gene, v_gene.split('*')[0]]:
        if vg in V_CONTRIBUTIONS:
            v_prefix = V_CONTRIBUTIONS[vg]
            break
    j_suffix = None
    for jg in [j_gene, j_gene.split('*')[0]]:
        if jg in J_CONTRIBUTIONS:
            j_suffix = J_CONTRIBUTIONS[jg]
            break
    if v_prefix is None or j_suffix is None:
        return None
    if not cdr3.startswith(v_prefix):
        for i in range(len(v_prefix), 1, -1):
            if cdr3.startswith(v_prefix[:i]):
                v_prefix = v_prefix[:i]
                break
        else:
            return None
    if not cdr3.endswith(j_suffix):
        for i in range(1, len(j_suffix)):
            if cdr3.endswith(j_suffix[i:]):
                j_suffix = j_suffix[i:]
                break
        else:
            return None
    n_start = len(v_prefix)
    n_end = len(cdr3) - len(j_suffix)
    if n_end < n_start:
        return None
    n_region = cdr3[n_start:n_end]
    return (v_prefix, n_region, j_suffix)

def get_bigrams(seq):
    return [seq[i:i+2] for i in range(len(seq)-1)]

# =============================================================================
# DATA LOADING
# =============================================================================

def load_vdjdb(chunks_dir='vdjdb-db/chunks'):
    records = []
    for fname in sorted(os.listdir(chunks_dir)):
        if not fname.endswith('.txt'):
            continue
        filepath = os.path.join(chunks_dir, fname)
        with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                records.append(row)
    return records

def filter_human_beta(records):
    return [r for r in records
            if r.get('cdr3.beta', '').strip()
            and r.get('species', '') == 'HomoSapiens'
            and r.get('v.beta', '').strip()
            and r.get('j.beta', '').strip()
            and r.get('antigen.epitope', '').strip()]

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

def analyze_epitope(records, epitope_seq, epitope_name=""):
    tcrs = [r for r in records if r.get('antigen.epitope') == epitope_seq]
    if not tcrs:
        print(f"No TCRs found for epitope {epitope_seq}")
        return {}
    print(f"
{'='*70}")
    print(f"EPITOPE: {epitope_seq} ({epitope_name})")
    print(f"{'='*70}")
    print(f"Total TCRs: {len(tcrs)}")
    decomposed = []
    for r in tcrs:
        result = decompose_cdr3(r['cdr3.beta'], r['v.beta'], r['j.beta'])
        if result:
            decomposed.append(result)
    success_rate = len(decomposed) / len(tcrs) * 100
    print(f"Successfully decomposed: {len(decomposed)} ({success_rate:.1f}%)")
    n_regions = [d[1] for d in decomposed if d[1]]
    n_counts = Counter(n_regions)
    print(f"
Top 15 N-region motifs (the 'semantic core'):")
    for nr, c in n_counts.most_common(15):
        pct = 100 * c / len(n_regions)
        bar = '█' * int(pct * 2)
        print(f"  '{nr}' — {c} times ({pct:.1f}%) {bar}")
    bigrams = Counter()
    for nr in n_regions:
        for bg in get_bigrams(nr):
            bigrams[bg] += 1
    total_bg = sum(bigrams.values())
    print(f"
Top 10 N-region bigrams ('phonotactics'):")
    for bg, c in bigrams.most_common(10):
        print(f"  {bg} — {100*c/total_bg:.2f}%")
    lengths = Counter(len(nr) for nr in n_regions)
    print(f"
N-region length distribution:")
    for l in sorted(lengths.keys())[:10]:
        bar = '█' * (lengths[l] // max(1, max(lengths.values()) // 30))
        print(f"  {l} aa: {lengths[l]:5d} {bar}")
    return {
        'epitope': epitope_seq, 'name': epitope_name,
        'total': len(tcrs), 'decomposed': len(decomposed),
        'n_regions': n_regions, 'n_counts': n_counts, 'bigrams': bigrams,
    }

def compare_epitopes(results_a, results_b):
    print(f"
{'='*70}")
    print(f"DISCRIMINATIVE ANALYSIS: {results_a['name']} vs {results_b['name']}")
    print(f"{'='*70}")
    bg_a, bg_b = results_a['bigrams'], results_b['bigrams']
    total_a, total_b = sum(bg_a.values()) or 1, sum(bg_b.values()) or 1
    ratios_a = {}
    for bg in bg_a:
        freq_a = bg_a[bg] / total_a
        freq_b = bg_b.get(bg, 0.5) / total_b
        ratios_a[bg] = freq_a / max(freq_b, 0.0001)
    print(f"
Bigrams enriched in {results_a['name']}:")
    for bg, ratio in sorted(ratios_a.items(), key=lambda x: -x[1])[:10]:
        print(f"  {bg} — {ratio:.1f}x enriched")
    ratios_b = {}
    for bg in bg_b:
        freq_b = bg_b[bg] / total_b
        freq_a = bg_a.get(bg, 0.5) / total_a
        ratios_b[bg] = freq_b / max(freq_a, 0.0001)
    print(f"
Bigrams enriched in {results_b['name']}:")
    for bg, ratio in sorted(ratios_b.items(), key=lambda x: -x[1])[:10]:
        print(f"  {bg} — {ratio:.1f}x enriched")

# =============================================================================
# MAIN
# =============================================================================

def main():
    print("Grammar of Immunity — Morphological Decomposition Demo")
    print("=" * 70)
    print("
Loading VDJdb...")
    records = load_vdjdb()
    print(f"Total records: {len(records)}")
    human_beta = filter_human_beta(records)
    print(f"Human TCRβ with full annotation: {len(human_beta)}")
    flu = analyze_epitope(human_beta, 'GILGFVFTL', 'Influenza A (M1)')
    cmv = analyze_epitope(human_beta, 'NLVPMVATV', 'CMV (pp65)')
    ebv = analyze_epitope(human_beta, 'GLCTLVAML', 'EBV (BMLF1)')
    if flu and cmv:
        compare_epitopes(flu, cmv)
    if flu and ebv:
        compare_epitopes(flu, ebv)
    print(f"
{'='*70}")
    print("DETAILED DECOMPOSITION EXAMPLES")
    print(f"{'='*70}")
    print()
    print(f"{'CDR3b':<25s} -> [{'V-prefix':<8s}] + [{'N-region':<10s}] + [{'J-suffix':<8s}]")
    print("-" * 70)
    examples = [
        ('CASSIRSSYEQYF', 'TRBV19', 'TRBJ2-7'),
        ('CASSIGAYGYTF', 'TRBV19', 'TRBJ1-2'),
        ('CASSGLAGLNEQFF', 'TRBV19', 'TRBJ2-1'),
        ('CASSPDQETSYTDTQYF', 'TRBV9', 'TRBJ2-3'),
        ('CSARLAGTEAFF', 'TRBV20-1', 'TRBJ1-1'),
    ]
    for cdr3, vg, jg in examples:
        result = decompose_cdr3(cdr3, vg, jg)
        if result:
            v, n, j = result
            print(f"{cdr3:<25s} -> [{v:<8s}] + [{n:<10s}] + [{j:<8s}]")
    print("

Done! See grammar-of-immunity-research.md for the full research plan.")

if __name__ == '__main__':
    main()
