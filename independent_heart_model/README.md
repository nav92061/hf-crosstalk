# Independent whole-heart model check (revision R1, response to R1.2.2)

Re-segmentation of the 87 primary-stratum series with the PlatiPy open whole-heart model
(Finnegan et al., Phys Eng Sci Med 2023; nnU-Net 3d_lowres, Zenodo record 6585664), trained on
radiotherapy heart contours from TCIA LCTSC and NSCLC-Radiomics rather than on TotalSegmentator output.

| File | Contents |
|---|---|
| `PRESPEC_independent_heart_model.md` | Analysis plan, written before any output existed |
| `run_platipy_heart.py` | GPU job: download each series from TCIA, segment, record checks and volume |
| `series.csv`, `series_pilot.csv`, `geometry.csv` | Series list, pilot pair, slice-position geometry used by the extent check |
| `volumes_platipy_heart.csv` | One row per series: heart volume, checks, timing, errors |
| `provenance_platipy.json` | Package versions, GPU, weight checksums, run times |
| `analyze_platipy.py` | Pre-specified analysis; reuses `imaging/endpoint1_tables.py` unchanged |
| `endpoint1_tables_ipp.py` | Same script with the slice-position z-extent (sensitivity analysis) |
| `independent_heart_model_results.csv` | Reported results |

Reproduce the analysis (CPU, seconds): `python independent_heart_model/analyze_platipy.py` from the repository root.
Re-running segmentation needs an NVIDIA GPU with platipy==0.7.2, nnunet==1.7.1, torch==2.0.1, vtk==9.3.1.
