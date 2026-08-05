# Pre-submission inquiry emails

A pre-submission inquiry asks an editor one question: *is this paper in scope for
your journal?* It is not a submission, it commits you to nothing, and you may send
several at once to different journals. The point is to learn in one week what a
blind submission would take six weeks to tell you — which matters here, because
the main risk to this paper is an editor deciding a negative result is not a
contribution.

## How to send them

- **Find the address.** Each journal page lists an editorial office email or an
  Editor-in-Chief. Springer journals usually have a "Contact the journal" link;
  Elsevier journals list the editorial office under "Editorial board" or in the
  guide for authors. If only a form exists, paste the text into the form.
- **Subject line:** `Pre-submission inquiry: TCR specificity benchmark, negative result`
- **Keep it under 250 words.** Editors skim. One paragraph on what you found, one
  on why it might fit, one asking the question.
- **Attach nothing**, or at most the abstract as plain text. Do not attach the
  manuscript.
- **Include the preprint DOI** once you have it. It lets an editor read the whole
  paper in one click if interested, which is the best possible outcome.
- **Send all four the same day.** Inquiries are not exclusive; submissions are. A
  positive reply obligates you to nothing, so there is no reason to stagger them
  or to hold any venue back.
- **Silence for two weeks means no.** Move on without following up.

## The question worth asking first

The single most decision-relevant unknown in the whole venue plan is whether
ImmunoInformatics will waive its USD 1,900 fee. Everything else about that journal
already favours it — it published both IMMREP workshop reports, so its editors and
reviewers are the community that built the benchmark you scored yourself against,
and a negative benchmarking result needs no justification there.

Cost is the only thing pushing it to fifth place in `paper/venue-plan.md`. If a
waiver is available, it goes to first. Ask on day one, not after three rejections.

## Reading the reply

- *"Please submit"* or *"this sounds within scope"* — submit there first.
- *"We would consider it but note that..."* — a real signal. Address the concern
  in the cover letter.
- *"Not within our scope"* — costs you nothing, saves you six weeks.
- Several rejections on **scope** grounds, rather than quality, is evidence the
  natural home is ImmunoInformatics and the USD 1,900 is buying a receptive
  audience rather than being a default cost.

---

## 1. Immunogenetics (Springer) — recommended first

> Dear Editors,
>
> We would like to ask whether the following study falls within the scope of
> *Immunogenetics*.
>
> We tested a published hypothesis: that segmenting TCR CDR3β along the boundaries
> created by V(D)J recombination — germline-contributed termini and the junctional
> N-region — yields receptor representations that generalize to epitopes absent
> from training. Using 121,467 deduplicated human TCRβ clonotypes from VDJdb and
> the official IMMREP23 benchmark scored with its official metric, we find that it
> does not. More informative for your readership, we find that the advantage such
> representations do show within the training distribution is explained almost
> entirely by germline V and J segment identity rather than junctional sequence,
> and that this holds across four negative-generation schemes, a six-point
> learning curve, and a retrained NetTCR-2.2.
>
> Two further results may interest immunogeneticists independently: a
> boundary-trim control showing that most apparent junctional structure sits at the
> germline boundaries, and a demonstration that a single tetramer-sort study
> contributing 79.5% of one epitope's clonotypes manufactures a convincing but
> artefactual specificity signature.
>
> All analyses are implemented in the Python standard library, are deterministic
> under fixed seeds, and are openly available. A preprint is openly deposited
> [DOI].
>
> Would this be of interest? We are unaffiliated researchers and would rather ask
> than presume.
>
> Sincerely,
> Erick Tejkowski (ORCID 0009-0006-9879-0777)
> Maria Elisa Paredes (ORCID 0009-0007-4967-8612)
> erick.tejkowski@gmail.com

## 2. Computers in Biology and Medicine (Elsevier)

Same study, framed as evaluation practice rather than immunogenetics.

> Dear Editors,
>
> We would like to ask whether the following study falls within the scope of
> *Computers in Biology and Medicine*.
>
> T-cell receptor specificity predictors perform well on epitopes seen during
> training and collapse to chance on unseen ones. We tested whether biologically
> informed receptor representations close that gap, on the official IMMREP23
> benchmark with its official metric. They do not: no method we evaluated exceeds
> 0.52 Macro AUC0.1 on the seven test peptides absent from training, including a
> published convolutional network retrained on the identical split, which reaches
> 0.92 AUC on its own validation pairs.
>
> The methodological findings may be the more transferable contribution. First, a
> six-point learning curve shows the representation's advantage rising
> monotonically with training size while a baseline stays flat — so a fixed-size
> comparison cannot distinguish a better representation from a more data-efficient
> one, which has implications for how representational claims in this field are
> evaluated. Second, we show how single-study composition in public databases can
> manufacture apparent biological signal.
>
> Every analysis is deterministic and reproducible from public code in minutes on a
> laptop. A preprint is openly deposited [DOI].
>
> Would a rigorous negative result of this kind be of interest?
>
> Sincerely,
> [as above]

## 3. Molecular Immunology (Elsevier)

> Dear Editors,
>
> We would like to ask whether the following study falls within the scope of
> *Molecular Immunology*.
>
> V(D)J recombination joins germline segments with non-templated junctional
> insertion, concentrating receptor diversity in CDR3. We asked whether the
> junctional region carries organization beyond composition, and whether exposing
> it improves prediction of which epitope a receptor recognizes. It does carry
> reproducible order-dependent structure across 79 held-out studies — though a
> boundary-trim control shows most of that signal sits at the germline boundaries
> rather than in the junctional interior. It does not improve prediction for
> epitopes absent from training, and the advantage seen for represented epitopes is
> explained by germline V and J identity rather than junctional sequence.
>
> We also document a data-quality issue with direct bearing on repertoire studies:
> for one of the three most-studied epitopes, a single study contributes 79.5% of
> clonotypes, and its junctional regions are 10.8% cysteine-containing against
> 0.4–0.6% elsewhere, producing a specificity signature that is experimental rather
> than immunological.
>
> All code and data are public and the analyses are deterministic. A preprint is
> openly deposited [DOI].
>
> Would this be of interest?
>
> Sincerely,
> [as above]

## 4. ImmunoInformatics (Elsevier) — send this one first

Send it with the others, on the same day. There is no reason to hold it back:
asking costs nothing, inquiries are not exclusive, and a "yes" obligates you to
nothing.

It is arguably the **most** valuable of the four, because it carries the one
question whose answer could reorder the entire plan. Scope is not in doubt — this
journal publishes the IMMREP post-mortems — so the real unknown is the fee. If a
waiver is available, ImmunoInformatics becomes both the best scientific fit and
free, and it goes from fifth choice to first. If no waiver is available, you have
learned that for the price of an email.

> Dear Editors,
>
> We have completed a direct test of a published hypothesis and believe
> *ImmunoInformatics* may be its natural home, given the journal's role in
> publishing the IMMREP benchmark reports.
>
> Vu et al. (*Nat Comput Sci*, 2024) argued that receptor language models fail
> partly because their tokenization encodes no immune-specific structure, and
> invited an implementation. We implemented it: CDR3β segmented into
> germline-contributed termini and the junctional N-region. On the official
> IMMREP23 data and metric, no method we evaluated exceeds 0.52 Macro AUC0.1 on the
> seven unseen peptides, including a retrained NetTCR-2.2, and the
> within-distribution advantage is explained by germline V/J identity rather than
> junctional sequence. Controls include a permuted-label capacity control, a
> size-matched negative control, a boundary-trim control, and a seen-epitope
> positive control at 0.71 AUC.
>
> Two questions. First, would the journal consider a negative result of this kind?
> Second, we are unaffiliated researchers with no funding of any kind, and we would
> like to ask whether an article publishing charge waiver may be available.
>
> A preprint is openly deposited [DOI].
>
> Sincerely,
> [as above]

---

## Before sending

Replace `[DOI]` in each with the preprint or Zenodo DOI once it arrives; see
`paper/preprint-deposit.md` for which comes first, and for why bioRxiv is no
longer the venue. Every number quoted above appears in the manuscript and
reproduces from the code — none of these emails claims anything the paper does
not support.
