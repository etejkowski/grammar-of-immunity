# Benchmark data and planned comparisons

## Included

`fetch.sh` downloads the official IMMREP23 challenge data from
https://github.com/justin-barton/IMMREP23 — training set (positives only),
test set, and `solutions.csv` with released labels. CSVs are gitignored.

Scoring uses the challenge's own metric, Macro AUC0.1 (per-peptide partial ROC
AUC to FPR 0.1, McClish-standardized, averaged). Verified against the challenge
documentation: an all-zero submission scores exactly 0.5000.

## Why NetTCR is not included, and what running it properly requires

NetTCR-2.2 (https://github.com/mnielLab/NetTCR-2.2) is the most tractable
published predictor to compare against: `src/predict.py` accepts exactly the
columns IMMREP23 supplies (`peptide`, `A1`, `A2`, `A3`, `B1`, `B2`, `B3`) and the
repository ships `.tflite` weights for fast inference.

It is nonetheless **not valid to score the IMMREP23 test set with those released
weights**:

* NetTCR-2.2's training data is compiled from IEDB, VDJdb and 10X Genomics.
* The repository additionally distributes IMMREP 2022 benchmark training data
  (`data/IMMREP/train/all_peptides.csv`).
* IMMREP23's TCRs are drawn from those same sources.

So the released model has likely seen a substantial share of the IMMREP23 test
TCRs during training. Any resulting score would be inflated by contamination —
and in the direction that flatters NetTCR against the models evaluated here,
which does not make it acceptable.

NetTCR-2.1's released models are peptide-specific. They cover six peptides
(ELAGIGILTV, GILGFVFTL, GLCTLVAML, IVTDFSVIK, NLVPMVATV, RAKFKQLL), of which
three occur in the IMMREP23 test set and **none** among its seven unseen
peptides. They therefore cannot speak to unseen-epitope generalization at all.

### The valid experiment

Retrain NetTCR-2.2 in `pan` mode on the identical IMMREP23 training split used
by `benchmark_immrep23.py`, with the same negative-generation protocol, then
score the official test set:

```
python src/train_nettcr_2_2_pan.py --train_data <our split> --val_data <held-out> ...
python src/predict.py --test_data <immrep23 test> --model_type pan ...
```

Requires TensorFlow (2.13–2.15 for Python 3.9). Expect a few hours including the
5-fold ensemble. This is the single most valuable addition to the manuscript and
is declared as future work in Limitations item 1.
