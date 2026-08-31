# Executable guardrails in two autonomous analysis pipelines

Reproduction material for the manuscript. Every number in the paper traces to a file here.

## Layout

    imaging/            imaging arm: cohort construction, headers, endpoints
    imaging/notebooks/  segmentation notebooks (Kaggle primary, Colab fallback)
    figures/            manuscript figures
    submission/         journal submission items

## Imaging arm, in order

1. `pull_headers.py`, one representative DICOM header per series from The Cancer Imaging
   Archive, via getSOPInstanceUIDs + getSingleImage. No volume transfer. Output:
   `headers_merged.csv` (573 series).
2. `pull_forid.py`, re-pull of FrameOfReferenceUID and AcquisitionTime for the 234 series
   in candidate pairs. This is what confirms two series share one irradiation event.
   Output: `recon_forid.csv`, `recon_pairs_confirmed.csv` (193 confirmed pairs).
3. `notebooks/phase1_segmentation_kaggle.ipynb`, TotalSegmentator 2.18.0 on the 217
   distinct series. Writes one row per series to `volumes.csv` and echoes each row to
   stdout prefixed `ROWCSV` so a lost output file is recoverable from the run log.
   Output: `volumes.csv`, `failure_rates.csv`, `model_provenance.json`.
4. `parse_run_log.py`, rebuilds the table from a pasted run log. Used for the reported
   run; the reconstructed table matched the recovered output file exactly.

## Endpoints

- Endpoint 1, reconstruction sensitivity: `endpoint1_clean_stratum.csv` is the primary
  result (45 coverage-matched pairs, 31 patients). `endpoint1_results.csv` carries the
  wider strata, `endpoint1_mixedmodel.csv` the second pre-specified aggregation.
- Endpoint 2, association stability: `endpoint2_verdict.json` records a null.
  `endpoint2_discordance_distribution.csv` integrates over the two arbitrary analyst
  choices; `endpoint2_null_discordance.csv` is the permutation calibration.
- Endpoint 3, failure rate: `endpoint3_failure_rate.json`. 4 of 217 series, all archive
  integrity faults arising before the model ran.

## Two things worth reading before the results

`endpoint1_prespec_violation.json` documents a pre-specification assumption that turned out
to be false. We had written that pairs sharing an acquisition therefore share z-coverage.
They do not: within-pair coverage differs by a median of 15%, and coverage difference
correlates with apparent volume difference at rho 0.70 for the aorta. The coverage-matching
filter that fixes this is post-hoc and is labelled as such.

`model_licensing_constraint.json` records why the structure list changed. The four cardiac
chambers originally pre-specified come from a task the TotalSegmentator package marks
commercial and gates behind a licence. The substitution was made on the open model before
any segmentation ran, and the amendment is dated in `PRESPECIFICATION_IMAGING.md`.

## Provenance

`model_provenance.json` pins TotalSegmentator 2.18.0, nnU-Net 2.8.1, PyTorch 2.10.0+cu128,
Tesla T4, and SHA-256 hashes for all six model checkpoints. Fast mode was not used: it
resamples to 3.0 mm, which would attenuate the thickness differences under study.

## Data

Nothing here is new data. Transcriptomic: GSE116250, GSE141910, and TCGA. Imaging:
TCGA-LUSC, TCGA-LUAD, TCGA-ESCA, TCGA-SARC in The Cancer Imaging Archive. Series and
instance identifiers are in `colab_series_manifest.csv`, so the exact scans are retrievable.
