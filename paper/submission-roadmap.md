# Step-by-step: from here to published

A chronological runbook. Every step says who does it, what to do, how long it
takes, and what to expect next. `paper/submission-packet.md` is the companion
document — it maps files to upload slots. This one is the order of operations.

Nothing below is urgent. The preprint step protects you regardless of how long
the rest takes.

---

## Stage 0 — Pre-flight (15 minutes, either author)

Confirm the deliverable still builds and passes its own checks.

```bash
cd /Users/e/grammar-of-immunity
python3 paper/build.py                     # DOCX + HTML + PDF + plain abstract
python3 paper/check_highlights.py          # expect PASS
.venv/bin/codespell -L FPR paper/*.md README.md   # expect no output
cp paper/build/manuscript.docx paper/TCR_Manuscript_V3.docx
```

Then open `paper/TCR_Manuscript_V3.docx` and read it once end to end. Not to
check the numbers — those are verified — but to hear the sentences. This is the
last point where changing your mind is free.

**Done when:** the DOCX opens, figures appear, and you have both read it.

---

## Stage 1 — Post the preprint (1 hour, Erick)

**bioRxiv is out.** It rejected the submission on 2026-08-04 (MS ID#
BIORXIV/2026/742614) because it requires authors to have an organizational
affiliation that can adjudicate ethical disputes. That is a policy gate applied
before anyone read the paper, and it changes nothing from Stage 4 onward.

The replacement route is in `paper/preprint-deposit.md`, which has the
field-by-field values for both forms:

1. **Zenodo–GitHub release** of the repository — no gatekeeper, code and data DOI
   the same day. `.zenodo.json` and `CITATION.cff` are written and waiting in the
   repo root; they must be committed and pushed before you cut the release,
   because Zenodo reads them out of the release tarball.
2. **OSF Preprints** deposit of `paper/build/manuscript.pdf` — the paper's
   preprint DOI. Free, indexed, no affiliation requirement, but moderated, and
   that file explains the one moderation risk to know about in advance.

Do them in that order and the inquiry emails never wait on a moderator.

**Done — 2026-08-05.** Both DOIs are live and resolving:

- preprint: `10.5281/zenodo.21800813` — https://zenodo.org/records/21800814
- code and data: `10.5281/zenodo.21800691` (release v1.0.0)

The preprint went to Zenodo rather than OSF Preprints; `paper/preprint-deposit.md`
explains why. Both DOIs are already in the manuscript, the cover letter, all four
inquiry letters, the README and `CITATION.cff`.

---

## Stage 2 — Screening (1–3 days, no action)

Zenodo does not screen; the DOI is immediate. OSF Preprints is pre-moderated,
usually a day or two, and tells you the reason if it declines. Neither is peer
review.

**What arrives:** a DOI, formatted like `10.5281/zenodo.xxxxxxx` from Zenodo or
`10.31219/osf.io/xxxxx` from OSF.

---

## Stage 3 — Record the DOI (5 minutes)

Four files need it:

1. `paper/cover-letter.md` — done, both DOIs in place
2. `paper/presubmission-inquiries.md` — done, all four letters carry both DOIs
3. `README.md` — done, badge and citation block at the top
4. `paper/manuscript.md` — done, "Data and code availability" now cites the
   Zenodo DOI alongside commit `53226f2`

Then `python3 paper/build.py` and re-copy the DOCX, so the manuscript and the
letters agree. The assistant can do this in one pass if you paste the DOI.

**Why it matters:** from this moment you have a timestamped, citable record with
both your names on it, independent of any journal's decision. Whatever happens
next, the work is yours and it is public.

---

## Stage 4 — Journal preparation (1 hour, Erick)

**Decide the venue first.** `paper/venue-plan.md` lays out the routes that cost
nothing, and recommends **Immunogenetics** (Springer) as the first submission:
free if you decline open access, indexed in PubMed, and its scope covers immune
receptor genetics, which is what the paper's central finding is about.
ImmunoInformatics remains the best community fit but is fully open access at USD
1,900 and is not in PubMed, so it is a deliberate paid choice rather than the
default.

Whichever you pick, open its guide for authors in a browser and find these:

1. Current publishing model, and whether any page or colour-figure charges apply
   even without open access. **Models flip** — two journals I had listed as free
   are no longer free.
2. Abstract word limit. Ours is 242 words, so a 250-word cap no longer bites.
   Only a cap of 200 or below forces a cut, starting with the TCRdist
   parenthetical.
3. Whether structured-abstract headings are fixed. Ours end with "Validation",
   which is non-standard.
4. Whether Highlights are required, and their character limit. Five are ready.
5. Reference style, and whether figure captions belong in the manuscript or a
   separate file.
6. Their generative-AI policy. Your declaration is already in the manuscript
   before the references; read how they frame it so you are not surprised later.

Paste those sections to the assistant; the files get conformed in one pass.

**Done when:** the venue is chosen, the abstract is within its limit, headings
match its format, and the highlights file is saved as .docx.

---

## Stage 5 — Submit to the journal (45 minutes, Erick)

Create an Editorial Manager account for ImmunoInformatics. Upload, in this order:

| Item | File |
|---|---|
| Manuscript | `paper/TCR_Manuscript_V3.docx` |
| Highlights | `paper/highlights.md` saved as .docx |
| Figures | `paper/figures/vector/*.pdf` |
| Cover letter | `paper/cover-letter.md`, adapted using the ImmunoInformatics note in that file |

In the cover letter, state the preprint DOI. Declare no competing interests and
no funding. If it asks for suggested reviewers, use the guidance at the bottom of
the cover letter, and do not name anyone either of you has corresponded with.

The system will generate a PDF proof of your submission. **Read it before
approving** — Editorial Manager mangles tables and special characters more often
than you would expect, and β and en-dashes are the usual casualties.

**Done when:** status reads "submitted to journal" or "with editor".

---

## Stage 6 — Waiting (2 weeks to 6 months)

Editorial Manager statuses, in the order you will probably see them:

| Status | Meaning |
|---|---|
| With Editor | An editor is deciding whether to send it out |
| Under Review | Reviewers have accepted the invitation |
| Required Reviews Complete | Reports are in; the editor is deciding |
| Decision in Process | A letter is coming |

DOAJ lists this journal's median submission-to-publication time at 23 weeks. A
desk rejection, if it happens, usually arrives within two weeks and means the
editor judged it out of scope rather than unsound.

**Do not** email to ask about progress before eight weeks. After twelve, a short
polite query is normal.

---

## Stage 7 — The decision

Four possible outcomes.

**Major revision.** The most likely, and it is a good outcome — it means they
engaged. Expect requests for one or more of: a comparison against TULIP or TCRen,
a VDJdb confidence-score replication, more detail on negative-example
construction, or a defense of the linear model class. All four are already
disclosed in Limitations, so none is a surprise. Go to Stage 8.

**Minor revision.** Clarifications and wording. Usually a week of work.

**Reject and resubmit.** Not the same as rejection. They want substantial new work
but would look again. Treat it as Stage 8 with a longer timeline.

**Rejection.** Go to Stage 9.

---

## Stage 8 — Revision (1–6 weeks depending on what they ask)

The mechanics that matter:

1. **Write a point-by-point response letter.** Quote each reviewer comment, then
   answer directly beneath it. Say exactly what changed and give the section and
   line. Reviewers look for whether you took them seriously.
2. **Where you disagree, say so with evidence, politely.** "We respectfully
   disagree, for the following reason" is normal and expected. Capitulating to a
   comment you believe is wrong makes the paper worse.
3. **Supply a tracked-changes version** as well as a clean one, if they ask.
4. **Run any new analysis in the repository**, not by hand, so the number in the
   response letter and the number in the code are the same thing.

The likely requests and their real cost, for planning:

| Request | Cost | Notes |
|---|---|---|
| VDJdb confidence-score replication | 1–2 days | needs the official release file, not the chunk records; will shift Table 1 slightly |
| Comparison against TULIP | 2–4 weeks | unsupervised, needs its own harness; no negative examples |
| Comparison against TCRen | 3+ weeks | needs structural models per pair |
| Clarify negatives / model class | 1 day | text only; the controls already exist |

---

## Stage 9 — If rejected

Nothing about a rejection invalidates the work, and your preprint stays up.

**Next venue, in my order of preference:** see `paper/venue-plan.md`, which plans
the route that costs nothing and corrects two journals I had wrongly listed as
free — Bioinformatics and Briefings in Bioinformatics both flipped to fully open
access, in 2023 and 2024 respectively. In short:

1. **Immunogenetics** (Springer) — free if you decline open access, PubMed-indexed,
   and its scope covers immune receptor genetics, which is what your central
   finding is about. This is where I would go first.
2. **Computers in Biology and Medicine** or **Molecular Immunology** (Elsevier) —
   also free without the OA option, also PubMed-indexed.
3. **PCI Mathematical and Computational Biology** → Peer Community Journal — free,
   rigorous, publishes the reviews alongside the paper, and requires the open code
   you already have. Weak PubMed presence is the tradeoff.
4. **Frontiers in Immunology** — widest readership, CHF 3,150, no fee support for
   US authors.

Whichever you choose, incorporate anything useful the first set of reviewers
said, even if they rejected it. Free peer review is still peer review.

---

## Stage 10 — Acceptance

1. **Settle any charges.** If you followed `paper/venue-plan.md` and submitted to
   a hybrid or subscription journal while declining open access, there is nothing
   to pay — decline the OA offer when it appears at acceptance, and check the
   acceptance letter for page or colour-figure charges. If you chose
   ImmunoInformatics, the APC is USD 1,900; ask about a waiver explicitly, since
   you are unaffiliated and unfunded.
2. **Proofs arrive within about two weeks.** Check the tables character by
   character — typesetters re-key numbers, and this paper is mostly numbers.
   Confirm β renders, the en-dashes in ranges survive, and every CI bracket is
   intact.
3. **Update the preprint.** Post the accepted version to whichever server holds
   the preprint as a new version; the original stays visible and the DOI still
   resolves. On Zenodo, a new GitHub release does the same for the code, and the
   concept DOI keeps resolving to the newest version.
4. **Add the published DOI** to the README and the repository description.

---

## Timeline, realistically

| Stage | Elapsed |
|---|---|
| Preprint live with a DOI | 1–3 days |
| Journal submission complete | within a week |
| First decision | 6 weeks to 6 months |
| Revision and resubmission | plus 2–8 weeks |
| Acceptance to publication | plus 4–8 weeks |

Roughly six to twelve months to a published paper, and about three days to a
citable one. That gap is the reason to do Stage 1 tonight and Stage 5 whenever
you are ready.
