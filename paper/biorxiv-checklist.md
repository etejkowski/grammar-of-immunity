# bioRxiv submission — what it asks for

A walkthrough so nothing is a surprise mid-form. Posting is free, takes roughly
20 minutes, and screening (a sanity check, not peer review) usually clears within
a couple of days.

Start at **https://www.biorxiv.org/submit-a-manuscript** and create an account
with the corresponding author's email.

## What to have ready

| Field | Our answer |
|---|---|
| Manuscript file | `paper/TCR_Manuscript_V3.docx` (figures embedded; built from `paper/manuscript.md`) |
| Figures | `paper/figures/*.png` — upload separately if prompted |
| Title | Biologically informed TCR representations improve in-distribution prediction but fail to generalize to unseen epitopes |
| Abstract | The Abstract section of the manuscript; paste as plain text |
| Corresponding author | Erick Tejkowski, erick.tejkowski@gmail.com |
| Authors, in order | Erick Tejkowski; Maria Elisa Paredes |
| ORCID | ET 0009-0006-9879-0777; MEP 0009-0007-4967-8612 |
| Institution | Leave blank or enter "Independent researcher" — bioRxiv does not require one |
| Subject category | **Bioinformatics** (alternative: Immunology) |
| Type | New Results |
| License | **CC-BY 4.0** recommended, or CC-BY-NC-ND if you want to bar commercial reuse |
| Competing interests | None |
| Funding | None |
| Data/code availability | https://github.com/etejkowski/grammar-of-immunity |
| Submitted to a journal? | No, at the point of posting |

## Things worth deciding before you start

**License.** CC-BY 4.0 lets anyone reuse with attribution and is what most
funders and journals prefer. It does not affect your ability to publish in a
journal afterward. CC-BY-NC-ND is more restrictive and occasionally creates
friction with journals later. I'd choose CC-BY 4.0.

**Author order.** ET first, MEP second, as in the manuscript. bioRxiv records
this permanently in the citation, so confirm it's what you both want.

**Both authors get an email.** bioRxiv notifies every listed author to confirm
authorship. Maria has to click hers or the posting stalls, so tell her to expect
it.

**Category choice.** Bioinformatics fits the methods and audience. Immunology
would reach immunologists but fewer methods-minded readers. Either is defensible;
Bioinformatics is my recommendation given the paper's emphasis on evaluation
practice.

## Sequence

1. Both ORCIDs are in hand and in the manuscript, so nothing blocks submission.
2. Post the preprint. You receive a DOI, typically within 1–2 days.
3. Add the DOI to `paper/cover-letter.md` where it says `DOI: XXXX`.
4. Submit to the journal, disclosing the preprint (all our target venues accept
   preprinted work).

## After posting

- The DOI is the citable, timestamped record — that's your protection as
  unaffiliated authors, independent of any journal timeline.
- bioRxiv pushes the preprint to your ORCID records automatically.
- Revisions post as new versions (v2, v3) with the original remaining visible.
- Add the DOI to the repository README so code and paper point at each other.

## What bioRxiv screening checks

Scope and plausibility, not quality: that it's research in the life sciences,
isn't plagiarized, contains no patient-identifying data, and makes no dangerous
claims. It is not peer review, and rejection at this stage is rare and usually
about category or format.

## What changed since the previous draft

- The NetTCR-2.2 retrain is **done**, not future work: retrained on the identical
  IMMREP23 split, it scores 0.6003 Macro AUC0.1 on seen peptides and 0.4868 on
  unseen, level with raw CDR3β 3-mers and below the TCRdist-style baseline. It is
  now a row in Table 5, a bar in Figure 3, and paragraphs in §2.7 and §3.6.
- In-text citations were added throughout; all 20 references are now cited.
- ET's ORCID was restored to the author block (it had been dropped) and MEP's
  real ORCID replaces "pending".
- `paper/manuscript.md` is again the single source of truth; the DOCX is built
  from it with `python3 paper/build.py`. Do not hand-edit the DOCX, or the two
  will diverge again.
