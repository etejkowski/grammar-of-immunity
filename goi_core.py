#!/usr/bin/env python3
"""
Grammar of Immunity — core data layer
=====================================

Shared by the phase scripts. Standard library only.

Responsibilities:

1. Load VDJdb chunk files.
2. Deduplicate clonotypes (VDJdb repeats the same (CDR3, V, J) across studies).
3. Derive germline V-prefix / J-suffix anchors *empirically* from the data,
   rather than from a hand-built table.
4. Decompose CDR3 into (V-prefix, N-region, J-suffix).

Why empirical anchors
---------------------
The original demo used a hand-typed table of ~40 V and ~13 J germline
contributions. That table was both incomplete (it silently dropped every
sequence using an unlisted segment, costing a third of the EBV set) and
wrong in a way that corrupted results (a strict endswith() test failed on
CDR3s carrying an extra trailing residue, so the J region leaked into the
N-region).

VDJdb bundles a partial reference at res/segments.aaparts.txt, but it covers
only 16 human TRBV and 13 TRBJ segments. Instead we derive the anchor for
every segment present in the data by consensus: for each position from the
5' end of CDR3s sharing a V gene, compute the modal residue and its
frequency; extend the prefix while modal frequency stays above a threshold.
The germline-encoded region is nearly invariant across clonotypes, while the
junctional region is not, so the frequency profile has a sharp drop at the
boundary. Same logic mirrored from the 3' end for J.

The reference file is then used to *check* the derived anchors, which is a
stronger position than trusting either source alone.
"""

import os
import csv
import math
from collections import Counter, defaultdict

CHUNKS_DIR = os.path.join('vdjdb-db', 'chunks')
AAPARTS = os.path.join('vdjdb-db', 'res', 'segments.aaparts.txt')

# Consensus threshold for calling a position "germline-encoded".
ANCHOR_THRESHOLD = 0.80
# Germline CDR3 contributions do not plausibly exceed this length.
MAX_ANCHOR = 8
# A segment needs at least this many distinct clonotypes to derive an anchor.
MIN_SEGMENT_SUPPORT = 20


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load_vdjdb(chunks_dir=CHUNKS_DIR):
    """Read every chunk TSV into a list of dicts."""
    records = []
    for fname in sorted(os.listdir(chunks_dir)):
        if not fname.endswith('.txt'):
            continue
        with open(os.path.join(chunks_dir, fname), encoding='utf-8',
                  errors='replace') as f:
            records.extend(csv.DictReader(f, delimiter='\t'))
    return records


def base_gene(g):
    """TRBV19*01 -> TRBV19."""
    return g.split('*')[0].strip()


VALID_AA = set('ACDEFGHIKLMNPQRSTVWY')


def filter_human_beta(records):
    """Human TCRbeta rows with the fields we need and clean sequences."""
    out = []
    for r in records:
        cdr3 = r.get('cdr3.beta', '').strip().upper()
        v = base_gene(r.get('v.beta', ''))
        j = base_gene(r.get('j.beta', ''))
        ep = r.get('antigen.epitope', '').strip().upper()
        if not (cdr3 and v and j and ep):
            continue
        if r.get('species', '') != 'HomoSapiens':
            continue
        if len(cdr3) < 6:
            continue
        if set(cdr3) - VALID_AA or set(ep) - VALID_AA:
            continue
        out.append({'cdr3': cdr3, 'v': v, 'j': j, 'epitope': ep,
                    'mhc': r.get('mhc.a', '').strip(),
                    'study': r.get('reference.id', '').strip() or 'unknown'})
    return out


def dedup_clonotypes(rows):
    """
    Collapse to unique (cdr3, v, j, epitope). VDJdb carries the same clonotype
    from multiple studies; counting rows inflates every downstream frequency.
    The row multiplicity is kept as 'n_rows' so it can be used as a weight,
    but analyses should default to one vote per clonotype. 'studies' retains
    the set of contributing studies, which is what makes cross-study
    replication testing possible.
    """
    merged = {}
    for r in rows:
        key = (r['cdr3'], r['v'], r['j'], r['epitope'])
        if key in merged:
            merged[key]['n_rows'] += 1
            merged[key]['studies'].add(r['study'])
        else:
            rec = dict(r)
            rec['n_rows'] = 1
            rec['studies'] = {r['study']}
            merged[key] = rec
    return list(merged.values())


# ---------------------------------------------------------------------------
# Empirical germline anchors
# ---------------------------------------------------------------------------

def derive_anchors(rows, threshold=ANCHOR_THRESHOLD, max_anchor=MAX_ANCHOR,
                   min_support=MIN_SEGMENT_SUPPORT):
    """
    Return (v_anchors, j_anchors, diagnostics).

    v_anchors[v_gene] = consensus 5' germline string
    j_anchors[j_gene] = consensus 3' germline string
    """
    by_v = defaultdict(list)
    by_j = defaultdict(list)
    for r in rows:
        by_v[r['v']].append(r['cdr3'])
        by_j[r['j']].append(r['cdr3'])

    def consensus(seqs, from_start):
        """Extend while the modal residue holds above threshold."""
        anchor = []
        profile = []
        for k in range(max_anchor):
            col = Counter()
            for s in seqs:
                if len(s) <= k + 1:      # never consume a whole CDR3
                    continue
                col[s[k] if from_start else s[-(k + 1)]] += 1
            if not col:
                break
            aa, n = col.most_common(1)[0]
            frac = n / sum(col.values())
            profile.append((aa, round(frac, 3)))
            if frac < threshold:
                break
            anchor.append(aa)
        seq = ''.join(anchor)
        return (seq if from_start else seq[::-1]), profile

    v_anchors, j_anchors, diag = {}, {}, {}
    for v, seqs in by_v.items():
        if len(seqs) < min_support:
            continue
        a, prof = consensus(seqs, True)
        if a:
            v_anchors[v] = a
            diag[('V', v)] = (len(seqs), a, prof)
    for j, seqs in by_j.items():
        if len(seqs) < min_support:
            continue
        a, prof = consensus(seqs, False)
        if a:
            j_anchors[j] = a
            diag[('J', j)] = (len(seqs), a, prof)
    return v_anchors, j_anchors, diag


def load_reference_anchors(path=AAPARTS, species='HomoSapiens', gene='TRB'):
    """
    Longest listed germline CDR3 contribution per segment from VDJdb's
    bundled reference. Partial coverage; used only to validate.
    """
    v_ref, j_ref = {}, {}
    if not os.path.exists(path):
        return v_ref, j_ref
    with open(path, encoding='utf-8', errors='replace') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['species'] != species or row['gene'] != gene:
                continue
            seg, part, typ = row['segm'], row['cdr3'], row['type']
            target = v_ref if typ == 'V' else j_ref
            if len(part) > len(target.get(seg, '')):
                target[seg] = part
    return v_ref, j_ref


# ---------------------------------------------------------------------------
# Decomposition
# ---------------------------------------------------------------------------

# Minimum informative anchor length. Unresolved gene names (bare families such
# as 'TRBV10' or 'TRBJ1') yield 2-residue consensus anchors like 'CA' or 'F',
# which carry no morphological information and would inflate the N-region.
MIN_ANCHOR_LEN = 3


def is_malformed(cdr3, j_anchors, min_len=4):
    """
    True for records whose 'CDR3' extends past the true CDR3 boundary.

    VDJdb contains a small number of entries where downstream framework (FR4)
    sequence was included, e.g. CSVPPGTDYNEQFFGPGTDYNEQFF, in which the
    germline J suffix appears twice. Such records break any morphological
    parse and must be excluded rather than parsed.
    """
    for s in j_anchors.values():
        if len(s) < min_len:
            continue
        first = cdr3.find(s)
        if first >= 0 and first + len(s) < len(cdr3) - 2:
            return True
    return False


def decompose(cdr3, v, j, v_anchors, j_anchors, min_j=2):
    """
    CDR3 -> (V-prefix, N-region, J-suffix) or None.

    The J anchor is located as the longest germline suffix occurring at or
    near the 3' end (tolerating up to 2 trailing residues), which is what the
    original strict endswith() test got wrong. Requiring at least `min_j`
    matched residues prevents a bare 'F' match from swallowing the J region
    into the N-region.
    """
    va = v_anchors.get(v)
    ja = j_anchors.get(j)
    if va is None or ja is None:
        return None
    if len(va) < MIN_ANCHOR_LEN or len(ja) < MIN_ANCHOR_LEN:
        return None

    # 5' side: longest prefix of the anchor that the CDR3 actually starts with.
    v_prefix = None
    for i in range(len(va), 0, -1):
        if cdr3.startswith(va[:i]):
            v_prefix = va[:i]
            break
    if not v_prefix:
        return None

    # 3' side: longest anchor suffix anchored near the end.
    j_suffix = None
    for i in range(0, len(ja) - min_j + 1):
        cand = ja[i:]
        idx = cdr3.rfind(cand)
        if idx >= 0 and 0 <= len(cdr3) - (idx + len(cand)) <= 2:
            j_suffix = cdr3[idx:]
            break
    if not j_suffix:
        return None

    n_start = len(v_prefix)
    n_end = len(cdr3) - len(j_suffix)
    if n_end < n_start:
        return None
    return v_prefix, cdr3[n_start:n_end], j_suffix


def annotate(rows, v_anchors, j_anchors):
    """
    Attach v_prefix / n_region / j_suffix; drop rows that will not parse.

    Returns (annotated_rows, stats) where stats records why rows were dropped.
    """
    out = []
    stats = {'malformed': 0, 'unparsed': 0}
    for r in rows:
        if is_malformed(r['cdr3'], j_anchors):
            stats['malformed'] += 1
            continue
        d = decompose(r['cdr3'], r['v'], r['j'], v_anchors, j_anchors)
        if not d:
            stats['unparsed'] += 1
            continue
        rec = dict(r)
        rec['v_prefix'], rec['n_region'], rec['j_suffix'] = d
        out.append(rec)
    return out, stats


def bigrams(seq):
    return [seq[i:i + 2] for i in range(len(seq) - 1)]


# ---------------------------------------------------------------------------
# Statistics (stdlib)
# ---------------------------------------------------------------------------

def _log_choose(n, k):
    if k < 0 or k > n:
        return float('-inf')
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1))


def fisher_exact_greater(a, b, c, d):
    """
    One-sided p-value for the 2x2 table [[a, b], [c, d]], testing whether a is
    over-represented. Pure hypergeometric tail sum.
    """
    row1, row2 = a + b, c + d
    col1 = a + c
    total = row1 + row2
    if total == 0 or row1 == 0 or col1 == 0:
        return 1.0
    logdenom = _log_choose(total, col1)
    p = 0.0
    kmax = min(row1, col1)
    for k in range(a, kmax + 1):
        lp = (_log_choose(row1, k) + _log_choose(row2, col1 - k) - logdenom)
        if lp > float('-inf'):
            p += math.exp(lp)
    return min(1.0, p)


def benjamini_hochberg(pvals):
    """Return BH-adjusted q-values, order preserved."""
    m = len(pvals)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [0.0] * m
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        i = m - rank + 1
        val = min(prev, pvals[idx] * m / i)
        q[idx] = val
        prev = val
    return q


def auc(scores_pos, scores_neg):
    """Mann-Whitney U / ROC AUC with tie handling."""
    labelled = [(s, 1) for s in scores_pos] + [(s, 0) for s in scores_neg]
    labelled.sort(key=lambda t: t[0])
    ranks = {}
    i = 0
    n = len(labelled)
    rank_sum_pos = 0.0
    while i < n:
        j = i
        while j + 1 < n and labelled[j + 1][0] == labelled[i][0]:
            j += 1
        avg_rank = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            if labelled[k][1] == 1:
                rank_sum_pos += avg_rank
        i = j + 1
    n_pos, n_neg = len(scores_pos), len(scores_neg)
    if n_pos == 0 or n_neg == 0:
        return float('nan')
    u = rank_sum_pos - n_pos * (n_pos + 1) / 2.0
    return u / (n_pos * n_neg)
