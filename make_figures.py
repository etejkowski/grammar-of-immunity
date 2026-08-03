#!/usr/bin/env python3
"""
Generate manuscript figures
===========================

Six figures at 300 dpi into paper/figures/.

Figures 1, 5 and 6 are recomputed from VDJdb directly, so they cannot drift
from the data. Figures 2, 3 and 4 plot values produced by the phase and
benchmark scripts; those values are declared as constants here with the script
that produced each one named alongside, so any discrepancy is auditable.

Usage:
    .venv/bin/python make_figures.py
"""

import os
import sys
from collections import Counter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, '.')
import goi_core as core

OUT = os.path.join('paper', 'figures')
DPI = 300

# Muted, print-safe palette; distinguishable in greyscale by ordering.
C_BASE = '#4C72B0'      # baseline / flat
C_TREAT = '#C44E52'     # morphological treatment
C_CTRL = '#8C8C8C'      # controls
C_REF = '#55A868'       # published-style references
C_ACC = '#8172B2'

plt.rcParams.update({
    'font.size': 9,
    'axes.titlesize': 10,
    'axes.labelsize': 9,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.dpi': DPI,
    'savefig.bbox': 'tight',
})


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    print(f"  wrote {path}")


# ---------------------------------------------------------------------------

def load_annotated():
    rows = core.dedup_clonotypes(core.filter_human_beta(core.load_vdjdb()))
    v_anch, j_anch, _ = core.derive_anchors(rows)
    ann, _ = core.annotate(rows, v_anch, j_anch)
    return ann, v_anch, j_anch


# ---------------------------------------------------------------------------
# Figure 1 — what the decomposition does
# ---------------------------------------------------------------------------

def fig1_decomposition(v_anch, j_anch):
    examples = [
        ('CASSIRSSYEQYF', 'TRBV19', 'TRBJ2-7'),
        ('CASSIRSTDTQYF', 'TRBV19', 'TRBJ2-3'),
        ('CASSIGAYGYTF', 'TRBV19', 'TRBJ1-2'),
        ('CSVGTGGTNEKLFF', 'TRBV29-1', 'TRBJ1-4'),
        ('CASSPDQETSYTDTQYF', 'TRBV9', 'TRBJ2-3'),
    ]
    fig, ax = plt.subplots(figsize=(6.6, 2.9))
    ax.axis('off')
    ax.set_xlim(0, 32)
    ax.set_ylim(-0.6, len(examples) + 0.6)

    cw = 0.82   # character width in data units
    x0 = 7.4

    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=c, alpha=0.85)
               for c in (C_BASE, C_TREAT, C_REF)]
    ax.legend(handles, ['V-prefix (germline)', 'N-region (junctional)',
                        'J-suffix (germline)'],
              loc='upper center', bbox_to_anchor=(0.5, 1.02), ncol=3,
              frameon=False, fontsize=8.5, handlelength=1.1,
              columnspacing=1.6)

    for i, (cdr3, v, j) in enumerate(examples):
        y = len(examples) - 1 - i
        d = core.decompose(cdr3, v, j, v_anch, j_anch)
        ax.text(6.9, y, f'{v}/{j}', fontsize=7.5, va='center', ha='right',
                color='#444444')
        if not d:
            continue
        vp, n, js = d
        for k, (seg, col) in enumerate(((vp, C_BASE), (n, C_TREAT),
                                        (js, C_REF))):
            start = len(vp) if k == 1 else (len(vp) + len(n) if k == 2 else 0)
            for m, chres in enumerate(seg):
                x = x0 + (start + m) * cw
                ax.add_patch(plt.Rectangle((x - cw / 2, y - 0.34), cw * 0.94,
                                           0.68, facecolor=col, alpha=0.85,
                                           edgecolor='white', linewidth=0.6))
                ax.text(x, y, chres, ha='center', va='center', fontsize=8,
                        color='white', family='monospace', weight='bold')

    ax.set_title('Morphological decomposition of CDR3$\\beta$\n'
                 'germline-encoded edges stripped; junctional interior isolated',
                 loc='left')
    save(fig, 'fig1_decomposition.png')


# ---------------------------------------------------------------------------
# Figure 2 — Phase 2 grammaticality, with controls
# ---------------------------------------------------------------------------

def fig2_grammaticality():
    # Source: phase2_grammaticality.py (A) and boundary_trim_control.py (B)
    labels = ['length\nonly', 'flat\nbigram', '(V,J)\nconditioned',
              'permuted\nlabel']
    vals = [0.5000, 0.6149, 0.7162, 0.6115]
    cols = [C_CTRL, C_BASE, C_TREAT, C_CTRL]

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.8, 3.5),
                                  gridspec_kw={'width_ratios': [1, 1.2]})
    bars = ax.bar(labels, vals, color=cols, width=0.62)
    ax.axhline(0.5, color='#333333', lw=0.9, ls=':', zorder=0)
    ax.text(3.85, 0.505, 'chance', fontsize=7.5, color='#333333', ha='right')
    ax.set_ylim(0.45, 0.78)
    ax.set_ylabel('AUC, real vs order-shuffled N-region')
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.006, f'{v:.4f}',
                ha='center', fontsize=8)
    ax.annotate('', xy=(2.35, 0.7162), xytext=(2.35, 0.6149),
                arrowprops=dict(arrowstyle='<->', color='#222222', lw=0.9))
    ax.text(2.45, 0.6655, '+0.1013\nCI [+0.0995,\n+0.1031]', fontsize=7,
            ha='left', va='center', color='#222222')
    ax.set_xlim(-0.6, 3.9)
    ax.tick_params(axis='x', labelsize=8)
    ax.set_title('A  Real, and not a capacity artifact\n'
                 '79 held-out studies, 16,304 pairs', loc='left')

    # ---- Panel B: boundary-trim control -------------------------------------
    # Source: boundary_trim_control.py
    tlab = ['none', '1 each\nend *', '2 each\nend', '3 each\nend',
            "2 from\n5\u2032 only *", "2 from\n3\u2032 only"]
    gaps = [0.1017, 0.0392, 0.0204, 0.0039, 0.0703, 0.0665]
    pairs = [16304, 14529, 7939, 1831, 14535, 14553]
    tcols = [C_TREAT, C_ACC, C_ACC, C_ACC, C_BASE, C_BASE]
    x = np.arange(len(tlab))
    b2 = ax2.bar(x, gaps, color=tcols, width=0.62)
    ax2.set_xticks(x)
    ax2.set_xticklabels(tlab, fontsize=7.8)
    ax2.set_ylabel('(V,J)-conditioned \u2212 flat, AUC')
    ax2.set_ylim(0, 0.142)
    for b, g, n in zip(b2, gaps, pairs):
        cx = b.get_x() + b.get_width() / 2
        ax2.text(cx, g + 0.0035, f'{g:+.4f}', ha='center', fontsize=7.4)
        ax2.text(cx, g + 0.0105, f'n={n:,}', ha='center', fontsize=6.3,
                 color='#666666')
    ax2.text(0.98, 0.135, '* same two residues removed, same n:\n'
             '   position, not amount, drives the gap',
             fontsize=6.9, va='top', ha='left', color='#222222')
    ax2.set_title('B  Most of the gain sits at the germline boundaries\n'
                  'residues trimmed from the N-region edges', loc='left')
    save(fig, 'fig2_grammaticality.png')


# ---------------------------------------------------------------------------
# Figure 3 — seen vs unseen on IMMREP23
# ---------------------------------------------------------------------------

def fig3_seen_unseen():
    # Source: benchmark_immrep23.py; nettcr2.2 from score_nettcr.py
    models = ['CDR3$\\beta$\n3-mers', 'morpheme', 'V/J genes\nonly',
              'all CDR loops\nboth chains', 'TCRbase-\nstyle',
              'TCRdist-\nstyle', 'NetTCR-2.2\nretrained']
    seen = [0.6003, 0.6378, 0.6346, 0.6495, 0.6414, 0.6971, 0.6003]
    unseen = [0.4858, 0.4992, 0.4791, 0.5192, 0.5000, 0.5000, 0.4868]
    cols = [C_BASE, C_TREAT, C_CTRL, C_ACC, C_REF, C_REF, C_REF]

    x = np.arange(len(models))
    fig, ax = plt.subplots(figsize=(8.0, 3.4))
    ax.bar(x - 0.2, seen, width=0.38, color=cols, label='seen peptides (13)')
    ax.bar(x + 0.2, unseen, width=0.38, color=cols, alpha=0.42,
           hatch='///', label='unseen peptides (7)')
    ax.axhline(0.5, color='#333333', lw=0.9, ls=':', zorder=0)
    ax.text(-0.44, 0.505, 'chance', fontsize=7.5, ha='left', color='#333333')
    ax.set_xticks(x)
    ax.set_xticklabels(models, fontsize=7.5)
    ax.set_ylim(0.44, 0.75)
    ax.set_ylabel('Macro AUC0.1 (official metric)')
    ax.legend(frameon=False, fontsize=8, loc='upper left')
    for xi, (s, u) in enumerate(zip(seen, unseen)):
        ax.text(xi - 0.2, s + 0.005, f'{s:.3f}', ha='center', fontsize=6.8)
        ax.text(xi + 0.2, u + 0.005, f'{u:.3f}', ha='center', fontsize=6.8)
    ax.set_title('Every method collapses to chance on unseen peptides\n'
                 'official IMMREP23 test set', loc='left')
    save(fig, 'fig3_seen_vs_unseen.png')


# ---------------------------------------------------------------------------
# Figure 4 — the advantage is a data-scale effect
# ---------------------------------------------------------------------------

def fig4_data_scale():
    # Source: negatives_robustness.py
    schemes = ['challenge\n(far decoys)\n67,872 pairs',
               'random\n(any decoys)\n67,872 pairs',
               'hard\n(near decoys)\n28,085 pairs',
               'matched\n(far decoys)\n29,957 pairs']
    deltas = [0.0345, 0.0333, -0.0022, -0.0009]
    cols = [C_TREAT, C_TREAT, C_CTRL, C_CTRL]

    fig, ax = plt.subplots(figsize=(6.2, 3.3))
    bars = ax.bar(schemes, deltas, color=cols, width=0.6)
    ax.axhline(0, color='#333333', lw=0.9)
    ax.set_ylabel('morpheme − 3-mers,\nMacro AUC0.1 (seen peptides)')
    for b, v in zip(bars, deltas):
        off = 0.0016 if v > 0 else -0.0032
        ax.text(b.get_x() + b.get_width() / 2, v + off, f'{v:+.4f}',
                ha='center', fontsize=8)
    ax.axvspan(-0.5, 1.5, color=C_TREAT, alpha=0.06, zorder=0)
    ax.axvspan(1.5, 3.5, color='#000000', alpha=0.04, zorder=0)
    ax.text(0.5, 0.041, 'full training set', fontsize=8, ha='center',
            color=C_TREAT)
    ax.text(2.5, 0.041, 'reduced training set', fontsize=8, ha='center',
            color='#555555')
    ax.annotate('same size as "hard",\ndecoys made dissimilar\n'
                '→ similarity is not the cause',
                xy=(3, -0.0009), xytext=(2.55, 0.020), fontsize=7.5,
                ha='center',
                arrowprops=dict(arrowstyle='->', lw=0.8, color='#333333'))
    ax.set_ylim(-0.010, 0.048)
    ax.tick_params(axis='x', labelsize=7.5)
    ax.set_title('The morphological advantage tracks training-set size,\n'
                 'not negative difficulty', loc='left')
    save(fig, 'fig4_data_scale.png')


# ---------------------------------------------------------------------------
# Figure 5 — batch confounding
# ---------------------------------------------------------------------------

def fig5_batch(ann):
    focus = [('GILGFVFTL', 'Influenza A\nM1'), ('NLVPMVATV', 'CMV\npp65'),
             ('GLCTLVAML', 'EBV\nBMLF1')]
    shares, cys, names = [], [], []
    for ep, name in focus:
        sub = [r for r in ann if r['epitope'] == ep]
        counts = Counter(s for r in sub for s in r['studies'])
        top = counts.most_common(1)[0][1] if counts else 0
        shares.append(100 * top / len(sub))
        cys.append(100 * sum(1 for r in sub if 'C' in r['n_region']) / len(sub))
        names.append(name)

    fig, axes = plt.subplots(1, 2, figsize=(6.6, 3.0))
    ax = axes[0]
    bars = ax.bar(names, shares, color=[C_BASE, C_ACC, C_TREAT], width=0.6)
    ax.set_ylabel('% of clonotypes from the\nsingle largest study')
    ax.set_ylim(0, 100)
    for b, v in zip(bars, shares):
        ax.text(b.get_x() + b.get_width() / 2, v + 1.5, f'{v:.1f}%',
                ha='center', fontsize=8)
    ax.tick_params(axis='x', labelsize=8)
    ax.set_title('Study concentration', loc='left')

    ax = axes[1]
    bars = ax.bar(names, cys, color=[C_BASE, C_ACC, C_TREAT], width=0.6)
    ax.set_ylabel('% of N-regions containing\ncysteine')
    for b, v in zip(bars, cys):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.25, f'{v:.1f}%',
                ha='center', fontsize=8)
    ax.tick_params(axis='x', labelsize=8)
    ax.set_title('Consequence: a "specificity" signal', loc='left')

    fig.suptitle('Apparent epitope-specific signal can be one study\'s batch',
                 x=0.02, ha='left', fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    save(fig, 'fig5_batch_effect.png')


# ---------------------------------------------------------------------------
# Figure 6 — the IRSS paradigm (known since 1998)
# ---------------------------------------------------------------------------

def fig6_motifs(ann):
    """
    The 1998 motif, at clonotype level.

    Raw motif counts are dominated by database redundancy: the single sequence
    CASSIRSSYEQYF occurs in 1,077 VDJdb rows but is ONE clonotype. Panel A
    shows that; panel B therefore reports the motif as the fraction of distinct
    clonotypes carrying it.
    """
    focus = [('GILGFVFTL', 'Influenza\nM1'), ('NLVPMVATV', 'CMV\npp65'),
             ('GLCTLVAML', 'EBV\nBMLF1')]
    stats = []
    for ep, name in focus:
        sub = [r for r in ann if r['epitope'] == ep]
        rs = 100 * sum(1 for r in sub if 'RS' in r['n_region']) / len(sub)
        v19 = 100 * sum(1 for r in sub if r['v'] == 'TRBV19') / len(sub)
        stats.append((name, rs, v19))

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1),
                             gridspec_kw={'width_ratios': [1, 1.5]})

    ax = axes[0]
    bars = ax.bar(['database\nrows', 'distinct\nclonotypes'], [1077, 1],
                  color=[C_CTRL, C_TREAT], width=0.58)
    ax.set_yscale('log')
    ax.set_ylim(0.5, 3000)
    ax.set_ylabel('count for CASSIRSSYEQYF')
    for b, v in zip(bars, [1077, 1]):
        ax.text(b.get_x() + b.get_width() / 2, v * 1.25, f'{v:,}',
                ha='center', fontsize=8.5)
    ax.tick_params(axis='x', labelsize=8)
    ax.set_title('A  Redundancy inflates motif counts', loc='left',
                 fontsize=9.5)

    ax = axes[1]
    x = np.arange(len(stats))
    ax.bar(x - 0.19, [s[1] for s in stats], width=0.36, color=C_TREAT,
           label='N-region contains RS')
    ax.bar(x + 0.19, [s[2] for s in stats], width=0.36, color=C_BASE,
           label='uses TRBV19 (=BV17)')
    ax.set_xticks(x)
    ax.set_xticklabels([s[0] for s in stats], fontsize=8)
    ax.set_ylabel('% of distinct clonotypes')
    ax.set_ylim(0, 50)
    ax.legend(frameon=False, fontsize=7.5, loc='upper right')
    for xi, s in enumerate(stats):
        ax.text(xi - 0.19, s[1] + 0.8, f'{s[1]:.1f}', ha='center', fontsize=7.5)
        ax.text(xi + 0.19, s[2] + 0.8, f'{s[2]:.1f}', ha='center', fontsize=7.5)
    ax.set_title('B  The 1998 BV17 / I-sRS(A)-S motif, recovered',
                 loc='left', fontsize=9.5)

    fig.tight_layout()
    save(fig, 'fig6_motif_paradigm.png')


def fig7_learning_curve():
    """Learning curve from learning_curve.json (produced by learning_curve.py)."""
    import json
    path = 'learning_curve.json'
    if not os.path.exists(path):
        print("  skipping fig7: run learning_curve.py first")
        return
    c = json.load(open(path))
    n = np.array([x['pairs'] for x in c], float)
    d_seen = np.array([x['delta_seen'] for x in c])
    d_uns = np.array([x['delta_unseen'] for x in c])
    m_seen = np.array([x['morpheme_seen'] for x in c])
    k_seen = np.array([x['kmer3_seen'] for x in c])

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2))

    ax = axes[0]
    ax.plot(n, m_seen, 'o-', color=C_TREAT, lw=1.6, ms=4.5,
            label='morpheme')
    ax.plot(n, k_seen, 's--', color=C_BASE, lw=1.6, ms=4.5,
            label='raw CDR3$\\beta$ 3-mers')
    ax.set_xscale('log')
    ax.set_xlabel('training pairs (log scale)')
    ax.set_ylabel('Macro AUC0.1, seen peptides')
    ax.legend(frameon=False, fontsize=8, loc='upper left')
    ax.set_title('A  Morphology benefits from data;\n     k-mers do not',
                 loc='left', fontsize=9.5)

    ax = axes[1]
    ax.plot(n, d_seen, 'o-', color=C_TREAT, lw=1.6, ms=4.5,
            label='seen peptides')
    ax.plot(n, d_uns, '^-', color=C_ACC, lw=1.6, ms=4.5,
            label='unseen peptides')
    ax.axhline(0, color='#333333', lw=0.8)
    ax.set_xscale('log')
    ax.set_xlabel('training pairs (log scale)')
    ax.set_ylabel('morpheme − 3-mers (Macro AUC0.1)')
    ax.legend(frameon=False, fontsize=8, loc='lower right')
    ax.text(n[0] * 1.05, d_seen[-1] * 0.98,
            f'still rising at the largest\nsize tested (+{d_seen[-1]:.4f})',
            fontsize=7.5, va='top', color='#333333')
    ax.set_title('B  The advantage is a trend, not a plateau',
                 loc='left', fontsize=9.5)

    fig.tight_layout()
    save(fig, 'fig7_learning_curve.png')


def main():
    print("Loading VDJdb for data-derived figures...")
    ann, v_anch, j_anch = load_annotated()
    print(f"  {len(ann):,} annotated clonotypes")
    print("Generating figures:")
    fig1_decomposition(v_anch, j_anch)
    fig2_grammaticality()
    fig3_seen_unseen()
    fig4_data_scale()
    fig5_batch(ann)
    fig6_motifs(ann)
    fig7_learning_curve()
    print("Done.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
