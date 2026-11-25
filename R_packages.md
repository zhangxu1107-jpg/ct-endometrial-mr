# R dependencies (code/MR_final_pipeline.R)

Required: `TwoSampleMR`, `dplyr`
Optional (the matching step is skipped if missing): `ggplot2` (forest plot), `MRPRESSO` (MR-PRESSO)

```r
install.packages(c("remotes", "dplyr", "ggplot2"))
remotes::install_github("MRCIEU/TwoSampleMR")
remotes::install_github("rondolab/MR-PRESSO")
```

Versions used for the reported results (from `sessionInfo()`):

- R 4.2.3 (2023-03-15), x86_64-apple-darwin17.0, macOS High Sierra 10.13.6
- TwoSampleMR 0.7.14, ieugwasr 1.2.1, dplyr 1.1.4, ggplot2 3.4.4
- MRPRESSO 1.0


