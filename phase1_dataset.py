#!/usr/bin/env python3
"""
Phase 1 — build a defensible dataset
====================================

Fixes the three methodological problems in the initial demo:

  1. hand-built germline table  -> empirical anchors, validated against
                                   VDJdb's bundled reference
  2. duplicate clonotype rows   -> deduplication on (CDR3, V, J, epitope)
  3. pseudocount enrichment     -> Fisher exact + Benjamini-Hochberg FDR,
                                   with a minimum-count floor

Outputs:
  annotated_cdr3.tsv   morphologically annotated, deduplicated dataset
  enrichment.tsv       bigram enrichment with q-values

Usage:
    python3 phase1_dataset.py
"""

import sys
from collections import Counter, defaultdict

import goi_core as core

MIN_BIGRAM_COUNT = 20      # count floor before a bigram is eligible
Q_THRESHOLD = 0.01         # FDR threshold for reporting

FOCUS = [
    ('GILGFVFTL', 'Influenza A (M1)'),
    ('NLVPMVATV', 'CMV (pp65)'),
    ('GLCTLVAML', 'EBV (BMLF1)'),
]


def rule(ch='='):
    print(ch * 74)


def main():
    rule()
    print("PHASE 1 — DATASET CONSTRUCTION")
    rule()

    raw = core.load_vdjdb()
    human = core.filter_human_beta(raw)
    uniq = core.dedup_clonotypes(human)
    print(f"VDJdb rows                        : {len(raw):>8,}")
    print(f"Human TCRbeta, fully annotated    : {len(human):>8,}")
    print(f"Unique clonotypes (cdr3,V,J,epi)  : {len(uniq):>8,}")
    dup_factor = len(human) / len(uniq) if uniq else 0
    print(f"Row inflation from duplication    : {dup_factor:>8.2f}x")

    # ---------------------------------------------------------------- anchors
    rule('-')
    print("EMPIRICAL GERMLINE ANCHORS")
    rule('-')
    v_anch, j_anch, _ = core.derive_anchors(uniq)
    print(f"V segments anchored: {len(v_anch)}   J segments anchored: {len(j_anch)}")

    v_ref, j_ref = core.load_reference_anchors()
    print(f"Reference coverage : {len(v_ref)} V, {len(j_ref)} J "
          f"(vdjdb-db/res/segments.aaparts.txt)")

    # Validate derived anchors against the reference where both exist.
    agree = disagree = 0
    mismatches = []
    for seg, derived in sorted(v_anch.items()):
        ref = v_ref.get(seg)
        if not ref:
            continue
        # The derived anchor should be a prefix of the reference contribution
        # (it may stop earlier, since we require 80% consensus).
        if ref.startswith(derived) or derived.startswith(ref):
            agree += 1
        else:
            disagree += 1
            mismatches.append(('V', seg, derived, ref))
    for seg, derived in sorted(j_anch.items()):
        ref = j_ref.get(seg)
        if not ref:
            continue
        if ref.endswith(derived) or derived.endswith(ref):
            agree += 1
        else:
            disagree += 1
            mismatches.append(('J', seg, derived, ref))
    print(f"Validation vs reference: {agree} consistent, {disagree} inconsistent")
    for typ, seg, derived, ref in mismatches[:10]:
        print(f"   MISMATCH {typ} {seg}: derived={derived!r} reference={ref!r}")

    print("\nSample anchors:")
    for seg in sorted(v_anch)[:6]:
        ref = v_ref.get(seg, '-')
        print(f"   V {seg:<10s} derived={v_anch[seg]:<8s} reference={ref}")
    for seg in sorted(j_anch)[:6]:
        ref = j_ref.get(seg, '-')
        print(f"   J {seg:<10s} derived={j_anch[seg]:<8s} reference={ref}")

    # ----------------------------------------------------------- decomposition
    rule('-')
    print("DECOMPOSITION")
    rule('-')
    ann, drop_stats = core.annotate(uniq, v_anch, j_anch)
    print(f"Decomposed: {len(ann):,} / {len(uniq):,} "
          f"({100 * len(ann) / len(uniq):.1f}%)")
    print(f"Dropped, malformed CDR3 (FR4 included): {drop_stats['malformed']:,}")
    print(f"Dropped, no usable anchor             : {drop_stats['unparsed']:,}")

    # Integrity assertion: no N-region may contain a full germline J suffix.
    j_full = {s for s in j_anch.values() if len(s) >= 4}
    leaks = [r for r in ann if any(s in r['n_region'] for s in j_full)]
    print(f"J-suffix leakage into N-region: {len(leaks)} "
          f"({'PASS' if not leaks else 'FAIL'})")
    if leaks:
        for r in leaks[:5]:
            print(f"   {r['cdr3']} {r['v']}/{r['j']} -> {r['n_region']}")
        return 1

    with open('annotated_cdr3.tsv', 'w', encoding='utf-8') as f:
        cols = ['cdr3', 'v', 'j', 'epitope', 'mhc', 'v_prefix', 'n_region',
                'j_suffix', 'n_rows']
        f.write('\t'.join(cols) + '\n')
        for r in ann:
            f.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')
    print("Wrote annotated_cdr3.tsv")

    # per-epitope summary
    print()
    print(f"{'epitope':<12s} {'source':<20s} {'rows':>7s} {'uniq':>7s} "
          f"{'parsed':>7s} {'rate':>7s}")
    for ep, name in FOCUS:
        rows_n = sum(1 for r in human if r['epitope'] == ep)
        uniq_n = sum(1 for r in uniq if r['epitope'] == ep)
        ann_n = sum(1 for r in ann if r['epitope'] == ep)
        rate = 100 * ann_n / uniq_n if uniq_n else 0
        print(f"{ep:<12s} {name:<20s} {rows_n:>7,} {uniq_n:>7,} "
              f"{ann_n:>7,} {rate:>6.1f}%")

    # ------------------------------------------------- study composition / batch
    rule('-')
    print("STUDY COMPOSITION (batch confounding check)")
    rule('-')
    print(f"{'epitope':<12s} {'clonotypes':>10s} {'studies':>8s} "
          f"{'largest study share':>20s}")
    for ep, name in FOCUS:
        sub = [r for r in ann if r['epitope'] == ep]
        counts = Counter(s for r in sub for s in r['studies'])
        top = counts.most_common(1)[0] if counts else ('-', 0)
        share = 100 * top[1] / len(sub) if sub else 0
        print(f"{ep:<12s} {len(sub):>10,} {len(counts):>8d} "
              f"{share:>18.1f}%  {top[0]}")
    print("\nAny epitope group dominated by a single study has its 'specific'"
          "\nsignal confounded with that study's batch. Enrichments below are"
          "\nreported with the number of independent studies replicating them.")

    # ------------------------------------------------------------- enrichment
    rule('-')
    print("BIGRAM ENRICHMENT — Fisher exact, BH-corrected, count floor "
          f"{MIN_BIGRAM_COUNT}")
    rule('-')

    def profile(ep, restrict_study=None):
        c = Counter()
        for r in ann:
            if r['epitope'] != ep:
                continue
            if restrict_study and restrict_study not in r['studies']:
                continue
            for bg in core.bigrams(r['n_region']):
                c[bg] += 1
        return c

    def studies_for(ep, min_clonotypes=200):
        counts = Counter(s for r in ann if r['epitope'] == ep
                         for s in r['studies'])
        return [s for s, n in counts.items() if n >= min_clonotypes]

    profiles = {ep: profile(ep) for ep, _ in FOCUS}
    out_rows = []
    for ep_a, name_a in FOCUS:
        for ep_b, name_b in FOCUS:
            if ep_a >= ep_b:
                continue
            for (x, nx), (y, ny) in ((((ep_a, name_a)), ((ep_b, name_b))),
                                     (((ep_b, name_b)), ((ep_a, name_a)))):
                pa, pb = profiles[x], profiles[y]
                ta, tb = sum(pa.values()), sum(pb.values())
                cand = [bg for bg, n in pa.items() if n >= MIN_BIGRAM_COUNT]
                pvals, recs = [], []
                for bg in cand:
                    a = pa[bg]
                    b = ta - a
                    c = pb.get(bg, 0)
                    d = tb - c
                    p = core.fisher_exact_greater(a, b, c, d)
                    ratio = ((a / ta) / ((c + 1) / (tb + 1)))
                    pvals.append(p)
                    recs.append((bg, a, c, ratio, p))
                qs = core.benjamini_hochberg(pvals)
                sig = [(r[0], r[1], r[2], r[3], r[4], q)
                       for r, q in zip(recs, qs) if q < Q_THRESHOLD]
                sig.sort(key=lambda t: -t[3])

                # cross-study replication for the reported bigrams
                x_studies = studies_for(x)
                per_study = {s: profile(x, s) for s in x_studies}
                per_study_tot = {s: sum(p.values()) or 1
                                 for s, p in per_study.items()}
                base_b = {bg: (pb.get(bg, 0) + 1) / (tb + 1) for bg in pa}

                print(f"\n{nx} > {ny}   ({len(sig)} significant of "
                      f"{len(cand)} tested; {len(x_studies)} studies with "
                      f">=200 clonotypes)")
                print(f"   {'bg':<4s} {'n_a':>6s} {'n_b':>6s} {'ratio':>8s} "
                      f"{'q':>10s} {'replicates':>11s}")
                for bg, a, c, ratio, p, q in sig[:8]:
                    rep = sum(1 for s in x_studies
                              if (per_study[s].get(bg, 0) / per_study_tot[s])
                              > base_b[bg])
                    flag = '' if not x_studies else (
                        f"{rep}/{len(x_studies)}")
                    print(f"   {bg:<4s} {a:>6d} {c:>6d} {ratio:>7.1f}x "
                          f"{q:>10.2e} {flag:>11s}")
                for bg, a, c, ratio, p, q in sig:
                    rep = sum(1 for s in x_studies
                              if (per_study[s].get(bg, 0) / per_study_tot[s])
                              > base_b[bg])
                    out_rows.append((nx, ny, bg, a, c, f"{ratio:.3f}",
                                     f"{p:.3e}", f"{q:.3e}",
                                     f"{rep}/{len(x_studies)}"))

    with open('enrichment.tsv', 'w', encoding='utf-8') as f:
        f.write("group_a\tgroup_b\tbigram\tcount_a\tcount_b\tratio\tp\tq\t"
                "studies_replicating\n")
        for row in out_rows:
            f.write('\t'.join(str(x) for x in row) + '\n')
    print("\nWrote enrichment.tsv")
    rule()
    return 0


if __name__ == '__main__':
    sys.exit(main())
