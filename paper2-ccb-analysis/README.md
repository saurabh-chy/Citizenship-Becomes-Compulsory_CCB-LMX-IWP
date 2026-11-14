# When Citizenship Becomes Compulsory — analysis code and outputs

This repository has the code and outputs for the manuscript *When Citizenship Becomes Compulsory: A Two-Wave Conditional-Process Study of Performance, Leader–Member Exchange, and Conscientiousness*.

The study uses two-wave survey data from 391 Indian IT/ITES employees, collected 8–9 months apart. It tests how compulsory citizenship behavior (CCB) relates to later leader–member exchange (LMX) and to individual work performance (IWP), and whether conscientiousness moderates these associations (H1–H7).

## Repository layout

```
code/
  paper2_analysis.ipynb   H1–H7, bootstrap mediation, CCB² sensitivity, Figures 1–7 and S4
  supporting_tables.py    descriptives, reliability, regression tables, diagnostics, VIFs, sensitivity checks
  measurement_lavaan.R    ordinal CFA, LMX invariance, CR/AVE, discriminant validity, latent SEM checks
outputs/
  paper2_main/            results_summary.txt, sensitivity_CCBsq.csv, figures/
  supporting_tables/      CSV tables written by supporting_tables.py
  measurement/            tables, figures, SEM outputs and session info from measurement_lavaan.R
data/                     not public; see data/README.md
```

## How to reproduce

Put the dataset in `data/` (see `data/README.md`), then run:

```bash
pip install -r requirements.txt
python code/supporting_tables.py --data data/CCB_Two_Phase_Data.csv --out outputs/supporting_tables
Rscript code/measurement_lavaan.R data/CCB_Two_Phase_Data.csv outputs/measurement
```

For the main analysis, open `code/paper2_analysis.ipynb` in Jupyter. By default it looks for `CCB_Two_Phase_Data.csv` in your Downloads folder; to use the copy in `data/`, edit the `DATA` line in the *Settings and data location* cell. Then run *Kernel → Restart & Run All*. Results and figures are written to a `Paper2_outputs` folder next to the data file.

The R script needs R ≥ 4.5 with `lavaan`, `semTools`, `semPlot`, `readr` and `psych`. Versions used are listed in `outputs/measurement/R_sessionInfo_complete.txt`.

The notebook takes under a minute. Its *Variables* cell checks N and three composite means against the values printed in the manuscript and flags any mismatch.

## Analysis settings

- **Scoring:** composites are item means, computed when at least 80% of a scale's items are answered. IWP = [task + contextual + (6 − CWB)] / 3.
- **Standardisation:** continuous predictors use z-scores with the population (N-denominator) SD. Outcomes stay on the 1–5 metric.
- **Controls:** gender, marital status, education (two dummies), total experience, organisational tenure and supervisor tenure.
- **Inference:** HC3 robust standard errors. Curvilinearity uses the Lind & Mehlum (2010) boundary-slope test.
- **Indirect effects:** computed as a × (b / SD of LMX2), with SD of LMX2 = 0.927 held fixed in every bootstrap and jackknife replicate. This puts both paths in raw LMX units.
- **Bootstrap:** 10,000 case resamples with seed 2026 and BCa intervals (jackknife acceleration). Each analysis starts its own generator from the same seed.

## Where each table and figure comes from

| Manuscript item | File(s) |
|---|---|
| Table 1 | `supporting_tables/Table_1_descriptives_correlations.csv` |
| Table 2 | `supporting_tables/Table_2_reliability.csv` |
| Table 3 | `paper2_main/results_summary.txt` (H3); `supporting_tables/Table_H3_hierarchical_strict_curvilinear_tests.csv` |
| Table 4 | `supporting_tables/Table_H6_quadratic_moderation.csv`, `Table_H6_conditional_curves.csv` |
| Table 5 | `paper2_main/sensitivity_CCBsq.csv` (primary rows) |
| Table 6 | `paper2_main/results_summary.txt` |
| Table S1 | `supporting_tables/Table_S3_participant_flow.csv`, `Table_S2_sample_descriptives.csv` |
| Table S2 | `supporting_tables/Table_1_descriptives_correlations.csv`, `Table_S9a_Pearson_correlations.csv`, `Table_S9b_Spearman_correlations.csv` |
| Table S3 | `supporting_tables/Table_2_reliability.csv`, `Table_S4_item_reliability_diagnostics.csv`; `measurement/tables/Table_R12_factorability_common_method.csv` |
| Table S4 | `measurement/tables/Table_R4_LMX_ordinal_invariance_FINAL.csv` |
| Table S5 | `supporting_tables/Table_4_HC3_quadratic_models.csv`, `Table_H1_H5_LMX_models.csv`, `Table_H2_LMX_IWP.csv` |
| Table S6 | `paper2_main/results_summary.txt` (H4, H5); `supporting_tables/Table_H5_conditional_slopes.csv` |
| Table S7 | `supporting_tables/Table_H6_quadratic_moderation.csv`, `Table_H6_conditional_curves.csv` |
| Table S8 | `paper2_main/sensitivity_CCBsq.csv` (H7 primary rows) |
| Table S9 | `measurement/tables/Table_R1_ordinal_CFA_fit.csv`, `Table_R13_CFA_estimator_sensitivity.csv` |
| Table S10 | `measurement/tables/Table_R9_CFA_CR_AVE.csv`, `Table_R21_complete_alpha_CR_AVE.csv` |
| Table S11 | `supporting_tables/Table_H3_hierarchical_strict_curvilinear_tests.csv` |
| Table S12 | `supporting_tables/Table_S11_control_sensitivity.csv`, `Table_S10_model_specific_VIF.csv` |
| Table S13 | `supporting_tables/Table_10_complete_diagnostics.csv` |
| Table S14 | `paper2_main/sensitivity_CCBsq.csv`, `results_summary.txt` (sensitivity block) |
| Figures 1–7, S4 | `paper2_main/figures/` (PNG and vector PDF) |

## Reproducibility notes

- **Main analysis:** `paper2_analysis.ipynb` produced the files in `outputs/paper2_main/` (original run: Python 3.8.8, matplotlib 3.3.4). Rerunning it on the study data reproduces `results_summary.txt` and `sensitivity_CCBsq.csv` exactly.
- **Figures:** the files in `paper2_main/figures/` are the ones in the manuscript. A rerun draws identical points, curves and confidence bands; on newer matplotlib versions, or machines without Times New Roman, the font and legend order can differ slightly.
- **Supporting tables:** `supporting_tables.py` reproduces all 21 CSVs in `outputs/supporting_tables/` exactly.
- **Measurement outputs:** the files in `outputs/measurement/` come from running `measurement_lavaan.R` on 4 September 2026 (seed 20260826).
