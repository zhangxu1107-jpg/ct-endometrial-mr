# Genetic proxies for *Chlamydia trachomatis* infection and endometrial cancer: Mendelian randomization code

Code and result tables for the manuscript "Genetic proxies for *Chlamydia trachomatis* infection and
endometrial cancer: a Mendelian randomization study with simulation-based assessment of instrument
validity" (submitted to *Infectious Disease Modelling*).

## Contents

| Path | Description |
|---|---|
| `code/MR_final_pipeline.R` | Full two-sample MR pipeline (reinfection and antibody exposures, sensitivity analyses, reverse MR, power, summary tables, forest plot) |
| `code/sim_duplicates.py` | Null simulation of duplicated instruments (400,000 replicates per scenario, seed 2026) |
| `data/Task1_exposure_Firth.csv` | The 22 suggestive reinfection SNPs (Firth beta, SE, P) taken from Tiwari et al. (2025) |
| `data/A_harmonised_13_small.csv` | Harmonised 13-SNP exposure-outcome data used by the simulation |
| `results/` | Result tables produced by the pipeline (`FINAL_summary.csv` holds every estimate reported in the paper) |

## Data sources

Summary statistics are not redistributed beyond the 22-SNP exposure table. They are obtained at run time:

- Reinfection exposure: Tiwari et al. (2025), African-ancestry cohort (283 women, 57 reinfected); SNPs in `data/`.
- Antibody exposure: OpenGWAS `ebi-a-GCST90006888` (*C. trachomatis* MOMP D antibody, n = 1,371).
- Endometrial cancer outcomes (OpenGWAS): `ebi-a-GCST90018838` (Sakaue et al., 2021; European part), `ebi-a-GCST006464`, `ebi-a-GCST006465`, `ebi-a-GCST006466` (O'Mara et al., 2018).

## How to run

1. Install the dependencies (`R_packages.md`, `requirements.txt`).
2. Get an OpenGWAS token at <https://api.opengwas.io> and set it **in your R session only**:
   ```r
   Sys.setenv(OPENGWAS_JWT = "<your token>")
   ```
   Do not write the token into any file in this repository (`.gitignore` excludes `.Renviron`).
3. From the repository root, run `source("code/MR_final_pipeline.R")`. Outputs are written to `final_out/`
   (an existing `final_out/` is renamed with a timestamp, not overwritten). The extraction test at the start
   stops the run with a clear message if the token is missing or expired.
4. Simulation: `python3 code/sim_duplicates.py data/A_harmonised_13_small.csv`
   (the pipeline also writes this input as `final_out/A_harmonised_13_small.csv`).

Expected simulation output (null, 400,000 replicates): 8 loci 3.8%; 13 SNPs 4.0%;
16 SNPs with 3 copies 8.9%; 13 SNPs plus 3 extra copies of the top-3 loci (22 rows) 23.7%.

## Notes

- OpenGWAS data and the API change over time; re-running may differ slightly from the reported tables if the underlying datasets are updated. The tables in `results/` are the outputs behind the manuscript.
- No credentials are stored in this repository.

## Citation

[Will insert the full citation and DOI of the paper after acceptance.]

## License

MIT (see `LICENSE`).
