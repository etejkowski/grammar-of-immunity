# Cover letter (journal submission)

Adapt the bracketed fields per venue. Keep it to one page — editors skim.

---

Dear Editors,

We submit for your consideration **"Biologically informed TCR representations
improve in-distribution prediction but fail to generalize to unseen epitopes."**

Jerne's 1984 Nobel lecture proposed that the immune system be understood as a
generative grammar, and Vu et al. (*Nature Computational Science*, 2024)
recently formalized that proposal for antibody sequences, arguing that current
receptor language models fail because their tokenization is linguistically
naive, and explicitly inviting an implementation. We implemented and tested it.

Our central finding is negative, and we believe it is useful precisely for that
reason. Morphological decomposition of CDR3β into germline-contributed edges and
the junctional interior does capture real structure: a (V,J)-conditioned model
separates real junctional regions from order-shuffled decoys at AUC 0.7162
versus 0.6149 for a flat model, on 79 entirely held-out studies, and a
permuted-label control confirms the gain is not model capacity. A boundary-trim
control we report alongside it shows most of that gain sits at the germline
boundaries, leaving +0.0392 attributable to the junctional interior — smaller,
but robust and, we think, the honest figure. But that
structure does not convert into predictive power for epitopes absent from
training. On the official IMMREP23 benchmark scored with the official Macro
AUC0.1 metric, no method we evaluated — including our own, reimplementations of
published approaches, and a retrained NetTCR-2.2 — exceeds 0.52 on the seven
test peptides absent from the training data.

Three further results may interest your readers beyond the immediate hypothesis:

1. The advantage of morphological tokenization over naive k-mers is **not
   separable from germline V/J segment usage** (+0.0091 Macro AUC0.1, CI
   crossing zero), corroborating existing concerns about V-gene shortcut
   learning by an independent route.
2. It is a **data-scale effect**: a six-point learning curve shows the advantage
   rising monotonically with training size and not plateauing, while naive
   k-mer performance stays flat. A fixed-size comparison cannot distinguish a
   better representation from a more data-efficient one, which we suggest has
   implications for how representational claims in this field are evaluated.
3. Apparent epitope-specific sequence signal in public data can be **dominated
   by single-study batch effects**: 79.5% of the GLCTLVAML clonotypes we
   analyzed originate from one study whose junctional regions are 10.8%
   cysteine-containing against 0.4–0.6% elsewhere, manufacturing a
   convincing-looking specificity signature.

To rule out the obvious objection that our null result reflects our model class,
we retrained NetTCR-2.2 — a published convolutional network over all six CDR
loops of both chains — on the identical IMMREP23 training split. It reaches
0.6003 Macro AUC0.1 on seen peptides and 0.4868 on unseen peptides, level with
raw CDR3β 3-mers (+0.0004, CI [−0.0425, +0.0339]) and significantly below a
TCRdist-style baseline (−0.0675, CI [−0.1164, −0.0272]), while fitting its own
validation pairs at AUC 0.9172. The generalization failure survives a change of
architecture.

We report our own corrected errors in the manuscript rather than omitting them,
including a row-versus-clonotype counting mistake and a single-seed false
positive that reseeding eliminated. The VDJdb analyses are implemented in the
Python standard library alone, are deterministic under fixed seeds, run in
minutes on a laptop, and are openly available at
https://github.com/etejkowski/grammar-of-immunity.

On our background, since we submit without institutional affiliation: ET's
formal training is in linguistics (MA) and software engineering, with immunology
grounding from pre-medical and dental study; MEP holds a BS in Biology, an MA in
Linguistics, and a PhD, and contributed the biological and immunological
interpretation. The linguistics-plus-immunology combination is the one Vu et al.
identified as necessary for this work and rarely available in one place.

This manuscript is not under consideration elsewhere. A preprint is openly
deposited at Zenodo, DOI 10.5281/zenodo.21800813, and the analysis code and data
pipeline are archived at DOI 10.5281/zenodo.21800691. We declare no competing interests and received no
funding.

Thank you for your consideration.

Sincerely,

Erick Tejkowski (ORCID 0009-0006-9879-0777)
Maria Elisa Paredes (ORCID 0009-0007-4967-8612)
Fairview Heights, Illinois, USA
erick.tejkowski@gmail.com

---

## Notes for adapting

- **ImmunoInformatics / Bioinformatics Advances** — lead with the benchmark result
  and the methodological cautions; both value reproducibility, so keep the
  stdlib-only and deterministic points.
- **PeerJ / BMC Bioinformatics** — these judge soundness rather than novelty, so
  emphasize the controls: positive control, capacity control, size-matched
  control, cross-study replication, seed stability, and the retrained NetTCR-2.2
  arm.
- **Frontiers in Immunology (Systems/Computational)** — lead with the immunology:
  the batch-effect finding and the germline-versus-junctional decomposition
  matter more to that readership than the tokenization question.
- If a venue requests **suggested reviewers**, plausible expertise sits with the
  IMMREP organizers, the ImmunoLingo group at Oslo, and authors of the 2025–2026
  TCR benchmarking assessments. Do not suggest anyone you have corresponded with.
- Some venues ask you to state why a negative result merits publication. The
  answer: the hypothesis was published and endorsed, the test was direct, the
  controls rule out the obvious confounds — including model class, now that a
  published deep architecture reproduces the failure — and the result redirects
  effort from the receptor side to the antigen side.
