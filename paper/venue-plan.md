# Where to publish without paying — venue plan

You do not pay reviewers. Reviewers are unpaid volunteers and editors are usually
academics working for free or a token stipend. The money I kept mentioning is an
**article publishing charge (APC)**: a fee some journals charge *authors* on
acceptance, to make the article free for everyone to read. It replaces
subscription income. There is no submission fee anywhere below, and nothing is
owed if a paper is rejected.

This document plans the route that costs nothing.

---

## Correction to earlier advice

I twice suggested Bioinformatics and Briefings in Bioinformatics as free options
because they were hybrid journals where you could decline open access. Both have
since flipped:

- **Bioinformatics** became fully open access in January 2023
- **Briefings in Bioinformatics** became fully open access in January 2024

Neither has a free route any more. Publishing models change, so **confirm the
model on the journal's own page before submitting** — that is exactly the mistake
I made, and it is easy to repeat.

## Why free publication still exists

Elsevier's own author-support page states it plainly: *"We are unable to offer
waivers for Article Publishing Charges in hybrid journals. These journals offer
the option to publish for free under a subscription model."*

So a **hybrid** journal gives you a choice: pay for open access, or publish at no
cost with the article behind a paywall. Since your preprint stays free on an open
server with its own DOI, anyone who wants the paper can read it either way. For this
project the paywall costs you very little.

The second free route is **diamond open access**: no charge to authors *or*
readers, funded by institutions rather than fees.

---

## Recommended order

**Revised 2026-08-05, after reading the editorial boards.** The earlier ranking
was built from cost, indexing and stated scope. Reading who actually runs each
journal, and what each has recently published, moved three of them. The old
ordering is preserved in the git history; what follows replaces it.

### 1. ImmunoInformatics (Elsevier) — best fit by a wide margin, USD 1,900 unless waived

Scope is no longer a question, and the fee is the only open one. Its 2025–2026
output is this paper's immediate neighbourhood:

- a tcrdist3 substitution-matrix comparison (now our reference [27])
- T-cell receptor specificity annotation models, June 2026
- active learning for out-of-distribution experimental design, March 2026
- "The gremlin in the works", on germline reference sequences, December 2025
  (now our reference [24])
- "Where single-cell transcriptomics fails T cells: the misuse of unsupervised
  clustering for T-cell annotation", December 2025

That last one is the decisive datum. It is a critical methodology paper, which
means the question "would you consider a negative result?" already has a
published answer. The germline paper is more consequential still: our manuscript
now engages it directly, and one of its authors, **Gur Yaari, is the journal's
Editor-in-Chief**.

The journal is also active, not winding down. Crossref records 62 articles: 16 in
2024, 17 in 2025, 3 by mid-2026. A library index listing coverage as 2021–2024
caused momentary alarm; that was the library's subscription window.

Not in PubMed, which remains a real cost to two unaffiliated authors building a
record. Ask about the waiver first; if granted, submit here.

### 2. Molecular Immunology (Elsevier) — free, PubMed-indexed, better fit than previously judged

**Promoted.** This document previously called the fit "looser" because the journal
leans molecular and cellular. That judgement predated knowing the editor.
**Zhinan Yin**, the Editor-in-Chief, is a T-cell immunologist whose recent work is
on Vγ9Vδ2 T-cell cytotoxicity, γδ T-cell tumour reactivity and TIL analysis. An
editor who publishes on TCR-bearing effector populations needs no persuading that
unseen-epitope generalization matters. 9,498 PubMed records.

### 3. Immunogenetics (Springer) — free, PubMed-indexed, fast

Unchanged on merit, and its scope still covers the finding that carries the paper:
that V and J identity, not junctional sequence, holds the transferable signal.
Two additions. Its median submission-to-first-decision is **9 days**, per
Springer's own metrics, which is exceptional and means little is lost by trying.
Springer publishes no editor email at all, so the contact came from a
corresponding-author footnote; see `paper/inquiry-emails/0-journal-contacts.txt`.

Hybrid, so declining open access should cost nothing. Verify before submitting.

### 4. Computers in Biology and Medicine (Elsevier) — free and well indexed, but a scope risk

**Demoted from second.** The cost and indexing arguments hold: free, 11,295 PubMed
records, and the evaluation-practice framing suits it. But its board has six
Executive Editors and no Editor-in-Chief, and not one of the six is an
immunologist or a sequence bioinformatician. The listed expertise is wavelets,
digital signal processing, biomedical image processing, computational
neuroscience, clinical decision support and drug repurposing. The closest match,
Sinosh Skariyachan, works in molecular modelling and computer-aided drug
discovery rather than repertoire analysis.

Desk rejection on scope is a live risk here, which is exactly why the inquiry has
value. Keep it in the first batch; do not count on it.

### 5. PCI Mathematical and Computational Biology → Peer Community Journal — free, and unusually well suited

`https://mcb.peercommunityin.org`

Unchanged. Free at every stage, diamond open access, reviews published alongside
the paper, and it requires the open code we already have — which suits work whose
main strength is its controls. The drawback is reach: 13 PubMed records, so
effectively unindexed. Check their policy on concurrent journal submission before
running it in parallel.

### 6. IEEE/ACM Transactions on Computational Biology and Bioinformatics — free, PubMed-indexed

Unchanged. 3,090 PubMed records, strong on methods and evaluation, but a more
computational audience than the paper's biology deserves. **Verify** whether page
charges are mandatory.

### 7. Frontiers in Immunology — CHF 3,150

Unchanged. Widest immunology readership, roughly USD 3,900, and their fee-support
programme excludes US authors. Only if reach justifies the cost.

---

## What reading the boards actually taught us

A pattern worth keeping in view when the replies arrive. The two immunology
journals are run by people whose own research is immune receptor biology; the
computational journal is not. If the replies split along that line — interest from
Immunogenetics and Molecular Immunology, scope doubts from Computers in Biology
and Medicine — that is evidence the paper's audience is immunologists who compute
rather than computational scientists who touch biology, and the eventual cover
letter should be framed accordingly.

One further connection: the ImmunoInformatics special issue carrying the germline
paper is guest-edited in part by **Justin Barton**, whose repository is the source
of the IMMREP23 data cited in our reference [12], and whose own nuTCRacker is now
our reference [26]. The people who curated the benchmark we scored ourselves
against are in that journal's editorial orbit. That favours us, but it means the
IMMREP23 methods will be read by someone who knows that data intimately.

---

## What I would actually do

**Send all four pre-submission inquiries on the same day** (drafts in
`paper/presubmission-inquiries.md`). They cost nothing, they are not exclusive, and
they answer in about a week what a blind submission takes six weeks to reveal.

Then let the replies decide:

- **If ImmunoInformatics grants a waiver**, submit there. Best fit, free, and the
  editor already published the concern our Limitations now engages.
- **If ImmunoInformatics says yes on scope but no on the fee**, the choice is
  USD 1,900 for the right audience against free publication in a PubMed-indexed
  immunology journal. Given that neither of us is funded, take the free route and
  keep the preprint doing the discoverability work.
- **Otherwise submit to Molecular Immunology or Immunogenetics**, whichever replies
  more warmly. Both are free, both are PubMed-indexed, and both are run by
  immunologists. Immunogenetics decides in a median of 9 days, so if both are
  positive, start there and lose almost nothing if it declines.
- **Treat Computers in Biology and Medicine as the fallback**, not the second
  choice. Its board has no immunologist on it.
- **If several decline on scope rather than quality**, that is evidence the natural
  audience is the IMMREP community, and paying the USD 1,900 becomes a considered
  choice rather than a default cost.

Run **PCI Mathematical and Computational Biology** in parallel or afterwards if you
want public, citable reviews at no cost. Check their policy on concurrent journal
submission first.

---

## Before submitting anywhere, confirm on the journal's own pages

1. Current publishing model — hybrid, fully OA, or subscription. Models flip.
2. Whether there are page, colour-figure or supplementary charges even without OA.
3. Abstract limit and heading format.
4. Whether highlights are required.
5. Their generative-AI policy, since the manuscript carries a declaration.

Paste those requirements to the assistant and every file is conformed in one pass.

## What does not change with venue

- The preprint is free and stays free, with a DOI and both your names.
- `paper/TCR_Manuscript_V3.docx`, the figures in PNG and vector PDF, the
  highlights, and the cover letter all work for any of these venues; only
  formatting details differ.
- The cover letter has per-venue notes at the bottom; the framing changes, the
  content does not.
