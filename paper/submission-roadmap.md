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
python3 paper/build.py                     # rebuilds DOCX + HTML from the markdown
python3 paper/check_highlights.py          # expect PASS
.venv/bin/codespell paper/*.md README.md   # expect no output
cp paper/build/manuscript.docx paper/TCR_Manuscript_V3.docx
```

Then open `paper/TCR_Manuscript_V3.docx` and read it once end to end. Not to
check the numbers — those are verified — but to hear the sentences. This is the
last point where changing your mind is free.

**Done when:** the DOCX opens, figures appear, and you have both read it.

---

## Stage 1 — Post the preprint (20 minutes, Erick)

Go to https://www.biorxiv.org/submit-a-manuscript and create an account with
erick.tejkowski@gmail.com.

Work through the form using the table in `paper/submission-packet.md`. The
answers that people get wrong:

- **Category:** Bioinformatics. Not Immunology. The paper's audience is
  methods-minded readers.
- **License:** CC-BY 4.0. It does not affect journal submission afterward.
- **Submitted to a journal?** No. That is true today and stays true until Stage 5.
- **Authors in order:** Erick, then Maria. bioRxiv freezes this into the citation.

Upload `paper/TCR_Manuscript_V3.docx`. If it asks for figures separately, use
`paper/figures/*.png`.

**Immediately after submitting, tell Maria to watch for an email from bioRxiv
asking her to confirm authorship.** The posting does not proceed until she clicks
it. This is the single most common cause of a stalled preprint.

**Done when:** you see "submitted" status and Maria has clicked her link.

---

## Stage 2 — Screening (1–2 days, no action)

A bioRxiv staff member checks that this is life-science research, is not
plagiarized, contains no patient-identifying data, and makes no dangerous claims.
It is not peer review. Rejection here is rare and is usually about category or
file format.

**What arrives:** an email with your DOI, formatted like
`10.1101/2026.08.03.xxxxxx`.

---

## Stage 3 — Record the DOI (5 minutes)

Two files need it:

1. `paper/cover-letter.md` — replace `DOI: XXXX`
2. `README.md` — add it near the top

Then `python3 paper/build.py` and re-copy the DOCX, so the manuscript and the
letter agree. The assistant can do this in one pass if you paste the DOI.

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
2. Abstract word limit. Ours is 264 words. If the cap is 250, cut the TCRdist
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
3. **Update the preprint.** Post the accepted version to bioRxiv as v2; the
   original stays visible and the DOI still resolves.
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
