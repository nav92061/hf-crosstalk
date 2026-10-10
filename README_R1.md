# Revision R1 additions (METRAD-D-26-00235)

Files added for the revised manuscript. Folder names match the paths cited in the paper.

| Folder | Contents | Reproduce |
|---|---|---|
| `imaging/` | `endpoint1_tables.py` (produces Table 2), its inputs, slice-position z-extent (`compute_z_extent.py`, `series_geometry_ipp.csv`, `slice_positions.csv.gz`), 47-pair table | `cd imaging && python endpoint1_tables.py` reproduces `endpoint1_clean_stratum.csv` exactly. `python compute_z_extent.py` re-downloads slice headers from TCIA (header data only, ~1 h). |
| `modules/` | The six analytical modules (`kernel.py` + `README.md` each). No language-model or platform calls. | Imported as plain Python files. |
| `plain_script/` | Agent-free reproduction of the pan-cancer ranking and guardrail tests (Table 5) | `python plain_script/run_plain.py plain_script/inputs modules/crosstalk_score/kernel.py out` (set `N_NULL=1000` for the published null) |
| `revision/` | Parameter audit, coverage sensitivity, agreement analysis, field-of-view sensitivity, CI and metric verification, chronology, threshold references | Result tables cited in Sections 2–3 |
| `agent_layer/` | Agent-layer record and all 23 task statements | Disclosure for Section 2.5 |
| `specifications/` | Pre-specification documents | |

Requirements: Python 3.10+, pandas, numpy, scipy, statsmodels, requests.
The agent runtime used to execute the analysis is a commercial application and is not included.
