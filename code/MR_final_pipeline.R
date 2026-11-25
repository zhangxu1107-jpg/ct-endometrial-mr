###############################################################################
# MR final pipeline: C. trachomatis (reinfection; MOMP D antibody) -> endometrial cancer
# Run once; all numbers reported in the paper come from the output of this script (FINAL_summary.csv)
#
# Requires: data/Task1_exposure_Firth.csv (reinfection exposure); internet access to OpenGWAS (token required)
# Packages: TwoSampleMR, dplyr (required); ggplot2, MRPRESSO (optional; the corresponding steps are skipped if missing)
#   MRPRESSO installation:remotes::install_github("rondolab/MR-PRESSO")
#   and that the antibody betas are in SD units (this affects R2).
###############################################################################
suppressMessages({ library(TwoSampleMR); library(dplyr) })
WD    <- ""                            # working directory (repository root); "" uses the current directory
TOKEN <- ""                            # OpenGWAS token: keep empty if set with Sys.setenv(OPENGWAS_JWT = ...); never commit a token
if (nzchar(WD)) { if (!dir.exists(WD)) stop("WD path does not exist: ", WD); setwd(WD) }
message("Working directory: ", getwd())
if (nzchar(TOKEN)) Sys.setenv(OPENGWAS_JWT = TOKEN)
if (nchar(Sys.getenv("OPENGWAS_JWT")) < 20) stop("No OpenGWAS token found: set TOKEN or run Sys.setenv(OPENGWAS_JWT = '...') first")
if (!file.exists("data/Task1_exposure_Firth.csv")) stop("Cannot find data/Task1_exposure_Firth.csv in: ", getwd())
options(warn = 1)
set.seed(2026)
# Do not overwrite or mix old results: archive an existing output folder, then create an empty final_out
if (dir.exists("final_out")) file.rename("final_out", paste0("final_out_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create("final_out")
OUT <- function(f) file.path("final_out", f)
cat("R version:", R.version.string, "| TwoSampleMR:", as.character(packageVersion("TwoSampleMR")), "\n")
invisible(tryCatch(ieugwasr::api_status(), error = function(e) NULL))   # silent call; the real token check is the data-extraction test below
.test <- tryCatch(extract_outcome_data("rs10208306", "ebi-a-GCST90018838", proxies = FALSE), error = function(e) NULL)
if (is.null(.test) || nrow(.test) == 0)
  stop("OpenGWAS extraction test failed: invalid or expired token, or a network problem. Generate a new token at https://api.opengwas.io, ",
       "set it with Sys.setenv(OPENGWAS_JWT = '<token>'), then run this script with source() (do not paste line by line).")
cat("OpenGWAS extraction test passed\n")
# Each part runs independently: an error in one part prints the reason and the script continues
step <- function(name, expr) tryCatch({ force(expr); cat("\n[", name, "] done\n") },
  error = function(e) cat("\n!!! [", name, "] error:", conditionMessage(e), "\n"))
outA <- list(); outB <- list()

## ---------------------------- Parameters ----------------------------------
ID_SMALL <- "ebi-a-GCST90018838"   # EC, smaller dataset (2,188 cases)
ID_LARGE <- "ebi-a-GCST006464"     # EC, larger dataset (12,906 cases)
ID_ENDO  <- "ebi-a-GCST006465"     # endometrioid subtype (8,758 cases; subset of the larger study)
ID_NONE  <- "ebi-a-GCST006466"     # non-endometrioid subtype (1,230 cases; subset of the larger study)
ID_AB    <- "ebi-a-GCST90006888"   # C. trachomatis MOMP D antibody
IDS    <- c(small = ID_SMALL, large = ID_LARGE, endo = ID_ENDO, nonendo = ID_NONE)
OUTLAB <- c(small = "EC small (2,188 cases)", large = "EC large (12,906 cases)",
            endo = "EC endometrioid (8,758)", nonendo = "EC non-endometrioid (1,230)")
N_EXPO <- 1371                      # sample size of the antibody GWAS
EXA <- "Ct reinfection (AFR, Tiwari)"
EXB <- "Ct MOMP D antibody (EUR)"
NCASE <- c(small = 2188, large = 12906, endo = 8758, nonendo = 1230)   # number of cases per outcome (from the dataset descriptions)

## ---------------------------- Helper functions ----------------------------
Isq <- function(b, se) { w <- 1/se^2; bm <- sum(w*b)/sum(w); Q <- sum(w*(b-bm)^2)
                         max(0, (Q-(length(b)-1))/Q) }
get_out <- function(snps, id, proxies) {
  o <- tryCatch(extract_outcome_data(snps, id, proxies = proxies), error = function(e) NULL)
  if (is.null(o) || nrow(o) == 0) NULL else o
}
harm <- function(expo, out, action) if (is.null(out)) NULL else harmonise_data(expo, out, action = action)
ok  <- function(d) if (is.null(d)) character(0) else d$SNP[d$mr_keep]
subs <- function(d, snps) if (is.null(d)) NULL else d[d$SNP %in% snps & d$mr_keep, , drop = FALSE]

RES <- list(); ZD <- list()
addres <- function(x) if (!is.null(x)) RES[[length(RES) + 1]] <<- x
est <- function(d, exposure, outcome, set, primary = FALSE, is_bin = TRUE) {
  if (is.null(d)) return(NULL)
  d <- d[d$mr_keep, , drop = FALSE]; n <- nrow(d); if (n == 0) return(NULL)
  if (n == 1) {
    b <- d$beta.outcome/d$beta.exposure; se <- d$se.outcome/abs(d$beta.exposure)
    p <- 2*pnorm(-abs(b/se)); Q <- NA; Qp <- NA; meth <- "Wald ratio"; i2 <- NA
  } else {
    r <- mr_ivw(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome)
    b <- r$b; se <- r$se; p <- r$pval; Q <- r$Q; Qp <- r$Q_pval; meth <- "IVW"
    i2 <- Isq(d$beta.exposure, d$se.exposure)
  }
  data.frame(exposure = exposure, outcome = unname(outcome), set = set, primary = primary,
    method = meth, nsnp = n, b = b, se = se, p = p,
    OR = if (is_bin) exp(b) else NA, OR_lo = if (is_bin) exp(b - 1.96*se) else NA,
    OR_hi = if (is_bin) exp(b + 1.96*se) else NA, Q = Q, Q_p = Qp, I2GX = i2,
    stringsAsFactors = FALSE)
}
addz <- function(label, r1, r2) {
  if (is.null(r1) || is.null(r2)) return(invisible())
  z <- (r1$b - r2$b)/sqrt(r1$se^2 + r2$se^2)
  ZD[[length(ZD) + 1]] <<- data.frame(comparison = label, z = z, p = 2*pnorm(-abs(z)))
}
wald_tab <- function(d) {
  if (is.null(d) || nrow(d) == 0)
    return(data.frame(SNP = character(0), b = numeric(0), se = numeric(0), p = numeric(0)))
  d <- d[d$mr_keep, ]
  b <- d$beta.outcome/d$beta.exposure; se <- d$se.outcome/abs(d$beta.exposure)
  data.frame(SNP = d$SNP, b = b, se = se, p = 2*pnorm(-abs(b/se)))
}
prune_loci <- function(expo_df, snps) {            # keep the SNP with the smallest Firth P within 500 kb
  cand <- expo_df %>% filter(rsID %in% snps) %>% arrange(FirthP); kept <- character(0)
  for (i in seq_len(nrow(cand))) { r <- cand[i, ]
    far <- vapply(kept, function(k) { kr <- cand[cand$rsID == k, ]
             !(kr$Chr == r$Chr && abs(kr$BP37 - r$BP37) < 500000) }, logical(1))
    if (all(far)) kept <- c(kept, r$rsID) }
  kept
}
power_tab <- function(R2s, ORs, nc, n, alpha = 0.05) {   # approximation for a binary outcome (Burgess 2014)
  nk <- n - nc; g <- expand.grid(R2 = R2s, OR = ORs)
  g$power <- mapply(function(R2, OR) { beta <- log(OR)
    se <- sqrt(1/R2 + beta^2)/sqrt(nc*nk/(nc + nk)); z <- abs(beta/se); za <- qnorm(1 - alpha/2)
    1 - pnorm(za - z) + pnorm(-za - z) }, g$R2, g$OR)
  g
}

###############################################################################
# A. Reinfection exposure (Tiwari, African ancestry): local table + online outcomes
###############################################################################
step("A", {
cat("\n################ A. Reinfection exposure ################\n")
expoA <- read.csv("data/Task1_exposure_Firth.csv", stringsAsFactors = FALSE)
cat("Firth P: <5e-8 =", sum(expoA$FirthP < 5e-8), "; <1e-6 =", sum(expoA$FirthP < 1e-6),
    "; <1e-5 =", sum(expoA$FirthP < 1e-5), "(of", nrow(expoA), "SNPs)\n")
cat("F (derived from P) range:", round(range(expoA$F_from_FirthP), 1), "\n")
fmtA <- function(df) { x <- suppressWarnings(suppressMessages(format_data(df, type = "exposure",
  snp_col = "rsID", beta_col = "beta_firth", se_col = "se_firth", pval_col = "FirthP",
  effect_allele_col = "EffectiveAllele", other_allele_col = "RefAllele")))
  x$exposure <- EXA; x$id.exposure <- "reinf"; x }

## A0. Unpruned 16-SNP analysis (includes 3 duplicate copies), smaller outcome, forward strand assumed
orig16 <- c("rs10208306","rs11996757","rs1669691","rs1669707","rs2486961","rs28505079","rs28530774",
            "rs66891172","rs6999003","rs7238797","rs73375993","rs7460431","rs7463208","rs7464228",
            "rs7829447","rs821932")
exp16 <- fmtA(expoA %>% filter(rsID %in% orig16))
d16   <- harm(exp16, get_out(exp16$SNP, ID_SMALL, FALSE), 1)
uniq13 <- expoA$rsID[expoA$rsID %in% orig16 & expoA$dup_keep == "keep"]
cat("Of the 16 SNPs, available in the smaller outcome:", length(ok(d16)), "\n")
write.csv(subs(d16, uniq13), OUT("A_harmonised_13_small.csv"), row.names = FALSE)   # input of sim_duplicates.py
addres(est(d16, EXA, OUTLAB["small"], "A0 unpruned analysis: 16 SNPs (with duplicate copies)"))
addres(est(subs(d16, uniq13), EXA, OUTLAB["small"], "A0 duplicate copies removed: 13 SNPs (forward strand assumed)"))
cat("[Reference] expected for the 16-SNP analysis: b about 0.0251, P about 0.0016\n")

## A1. Main flow: remove exact duplicates (18 SNPs) -> extract outcomes -> locus pruning
expo18 <- expoA %>% filter(dup_keep == "keep"); exp18 <- fmtA(expo18)
outA <- lapply(IDS, function(id) get_out(exp18$SNP, id, FALSE))     # no proxy SNPs (African-ancestry exposure)
dA3 <- lapply(outA, function(o) harm(exp18, o, 3))                  # exclude palindromic SNPs
dA1 <- lapply(outA, function(o) harm(exp18, o, 1))                  # same strand assumed, palindromic SNPs kept
avail3 <- intersect(ok(dA3$small), ok(dA3$large)); kept3 <- prune_loci(expo18, avail3)
avail1 <- intersect(ok(dA1$small), ok(dA1$large)); kept1 <- prune_loci(expo18, avail1)
cat("\nFlow: 22 -> ", nrow(expo18), " (exact duplicates removed) -> available in both outcomes: palindromic excluded ", length(avail3),
    " / palindromic kept ", length(avail1), " -> separate loci: ", length(kept3), " / ", length(kept1), "\n", sep = "")
cat("Separate loci (palindromic excluded):", kept3, "\nSeparate loci (palindromic kept):", kept1, "\n")

s3 <- l3 <- s1 <- l1 <- NULL
for (o in names(IDS)) {
  r3 <- est(subs(dA3[[o]], kept3), EXA, OUTLAB[o], "A separate loci (palindromic SNPs excluded)", primary = TRUE)
  r1 <- est(subs(dA1[[o]], kept1), EXA, OUTLAB[o], "A separate loci (palindromic SNPs kept)")
  addres(r3); addres(r1)
  if (o == "small") { s3 <- r3; s1 <- r1 }; if (o == "large") { l3 <- r3; l1 <- r1 }
}
for (o in c("small", "large"))
  addres(est(subs(dA3[[o]], avail3), EXA, OUTLAB[o], "A sensitivity: 18 de-duplicated SNPs (not pruned)"))
addz("A separate loci (palindromic excluded): small vs large", s3, l3); addz("A separate loci (palindromic kept): small vs large", s1, l1)
write.csv(wald_tab(subs(dA3$small, kept3)), OUT("A_wald_small_loci.csv"), row.names = FALSE)
write.csv(wald_tab(subs(dA3$large, kept3)), OUT("A_wald_large_loci.csv"), row.names = FALSE)

## A. Table S1
pal <- c(dA1$small$SNP[dA1$small$palindromic], dA1$large$SNP[dA1$large$palindromic])
tabA <- expoA %>% mutate(queried = dup_keep == "keep",
  in_small = rsID %in% outA$small$SNP, in_large = rsID %in% outA$large$SNP,
  palindromic = rsID %in% pal, in_loci_excl_pal = rsID %in% kept3, in_loci_keep_pal = rsID %in% kept1)
write.csv(tabA, OUT("TableS1_A_reinfection_SNPs.csv"), row.names = FALSE)
for (nm in names(IDS)) { if (!is.null(outA[[nm]])) write.csv(outA[[nm]], OUT(paste0("A_out_", nm, ".csv")), row.names = FALSE) }

})

###############################################################################
# B. Antibody exposure (MOMP D, European ancestry): online instruments + online outcomes
###############################################################################
step("B", {
cat("\n################ B. Antibody exposure ################\n")
snp13 <- c("rs117664663","rs2739625","rs556141332","rs8053331","rs186551975","rs11198255",
           "rs11571815","rs35775645","rs139532170","rs144941477","rs75572600","rs75732672","rs28694248")
xB_all <- extract_instruments(ID_AB, p1 = 5e-8, clump = FALSE)
xB13 <- xB_all[xB_all$SNP %in% snp13, ]
if (nrow(xB13) != 13) warning("Expected 13 instruments but got ", nrow(xB13), "; check before continuing")
instB <- xB13 %>% transmute(SNP, eaf = eaf.exposure, beta = beta.exposure, se = se.exposure,
  F = (beta/se)^2, R2 = 2*eaf*(1 - eaf)*beta^2, minor_alleles = round(2*N_EXPO*pmin(eaf, 1 - eaf)))
cat("Mean F =", round(mean(instB$F), 2), "; minimum F =", round(min(instB$F), 1),
    "; sum of R2 (13 SNPs) =", round(sum(instB$R2), 3), "(inflated by winner's curse; betas must be in SD units)\n")

outB <- lapply(IDS, function(id) get_out(xB13$SNP, id, TRUE))       # default proxy search kept
dB   <- lapply(outB, function(o) harm(xB13, o, 2))
s13 <- ok(dB$small); common <- intersect(s13, ok(dB$large)); rare3 <- setdiff(s13, ok(dB$large))
cat("Available in the smaller outcome:", length(s13), "; shared by both outcomes:", length(common), "; only in the smaller outcome:", rare3, "\n")
cat("Sum of R2 of the shared SNPs =", round(sum(instB$R2[instB$SNP %in% common]), 3), "\n")

rB <- list()
rB$s13  <- est(dB$small, EXB, OUTLAB["small"], "B all available SNPs", primary = TRUE)
rB$sc   <- est(subs(dB$small, common), EXB, OUTLAB["small"], "B SNPs shared by both outcomes", primary = TRUE)
rB$lc   <- est(subs(dB$large, common), EXB, OUTLAB["large"], "B SNPs shared by both outcomes", primary = TRUE)
rB$ec   <- est(subs(dB$endo, common), EXB, OUTLAB["endo"], "B SNPs shared by both outcomes (subtype, available SNPs)", primary = TRUE)
rB$nc   <- est(subs(dB$nonendo, common), EXB, OUTLAB["nonendo"], "B SNPs shared by both outcomes (subtype, available SNPs)", primary = TRUE)
rB$r3   <- est(subs(dB$small, rare3), EXB, OUTLAB["small"], "B rare variants available only in the smaller outcome")
rB$no1  <- est(subs(dB$small, setdiff(s13, "rs117664663")), EXB, OUTLAB["small"], "B all available SNPs without rs117664663")
for (r in rB) addres(r)
addz("B shared SNPs: small vs large", rB$sc, rB$lc)

sing <- mr_singlesnp(dB$small[dB$small$mr_keep, ]); loo <- mr_leaveoneout(dB$small[dB$small$mr_keep, ])
write.csv(sing, OUT("B_singleSNP_small.csv"), row.names = FALSE)
write.csv(loo,  OUT("B_leaveoneout_small.csv"), row.names = FALSE)

## B. Heterogeneity / pleiotropy / Egger / weighted median
tryCatch({
  for (nm in c("small", "large")) {
    dd <- if (nm == "small") subs(dB$small, s13) else subs(dB$large, common)
    cat("\n--", nm, "heterogeneity / Egger intercept --\n"); print(mr_heterogeneity(dd)); print(mr_pleiotropy_test(dd))
    print(mr(dd, method_list = c("mr_egger_regression", "mr_weighted_median"))) }
}, error = function(e) cat("Heterogeneity/pleiotropy step skipped: ", conditionMessage(e), "\n"))

## B. MR-PRESSO (larger outcome, shared SNPs)
tryCatch({
  dl <- subs(dB$large, common)
  pr <- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
          SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE,
          DISTORTIONtest = TRUE, data = dl, NbDistribution = 10000, SignifThreshold = 0.05)
  print(pr$`Main MR results`); cat("Global test P =", pr$`MR-PRESSO results`$`Global Test`$Pvalue, "\n")
  capture.output(print(pr), file = OUT("B_MRPRESSO_large.txt"))
}, error = function(e) cat("MR-PRESSO skipped (package not installed or error): ", conditionMessage(e), "\n"))

## B. Forward Steiger test (antibody -> EC, larger outcome, shared SNPs)
tryCatch({
  ds <- subs(dB$large, common); ds$samplesize.exposure <- N_EXPO
  print(directionality_test(ds))
  st <- steiger_filtering(ds); print(st[, c("SNP", "rsq.exposure", "rsq.outcome", "steiger_dir", "steiger_pval")])
  write.csv(st, OUT("B_steiger_forward.csv"), row.names = FALSE)
}, error = function(e) cat("Forward Steiger skipped: ", conditionMessage(e), "\n"))

## B. Sensitivity 1: no proxy SNPs
tryCatch({
  oS <- get_out(xB13$SNP, ID_SMALL, FALSE); oL <- get_out(xB13$SNP, ID_LARGE, FALSE)
  dS0 <- harm(xB13, oS, 2); dL0 <- harm(xB13, oL, 2); c0 <- intersect(ok(dS0), ok(dL0))
  addres(est(dS0, EXB, OUTLAB["small"], "B sensitivity: no proxies, all available SNPs"))
  addres(est(subs(dS0, c0), EXB, OUTLAB["small"], "B sensitivity: no proxies, shared SNPs"))
  addres(est(subs(dL0, c0), EXB, OUTLAB["large"], "B sensitivity: no proxies, shared SNPs"))
}, error = function(e) cat("No-proxy sensitivity analysis skipped: ", conditionMessage(e), "\n"))

## B. Sensitivity 2: instrument set returned by online clumping in extract_instruments() (11 SNPs)
tryCatch({
  xB11 <- extract_instruments(ID_AB, p1 = 5e-8, clump = TRUE, r2 = 0.001, kb = 10000)
  cat("Online clumping returned", nrow(xB11), "SNPs; missing from the 13:", setdiff(snp13, xB11$SNP), "; added:", setdiff(xB11$SNP, snp13), "\n")
  d11s <- harm(xB11, get_out(xB11$SNP, ID_SMALL, TRUE), 2); d11l <- harm(xB11, get_out(xB11$SNP, ID_LARGE, TRUE), 2)
  c11 <- intersect(ok(d11s), ok(d11l))
  addres(est(d11s, EXB, OUTLAB["small"], "B sensitivity: online clumping (11 SNPs), all available"))
  addres(est(subs(d11s, c11), EXB, OUTLAB["small"], "B sensitivity: online clumping (11 SNPs), shared SNPs"))
  addres(est(subs(d11l, c11), EXB, OUTLAB["large"], "B sensitivity: online clumping (11 SNPs), shared SNPs"))
}, error = function(e) cat("Online-clumping sensitivity analysis skipped: ", conditionMessage(e), "\n"))

## B. Table S1
tabB <- instB %>% mutate(
  in_small = SNP %in% outB$small$SNP, in_large = SNP %in% outB$large$SNP,
  proxy_small = SNP %in% outB$small$SNP[which(outB$small$proxy.outcome %in% TRUE)],
  proxy_large = SNP %in% outB$large$SNP[which(outB$large$proxy.outcome %in% TRUE)],
  only_small = SNP %in% rare3, palindromic = SNP %in% c(dB$small$SNP[dB$small$palindromic %in% TRUE],
                                                       dB$large$SNP[dB$large$palindromic %in% TRUE]))
write.csv(tabB, OUT("TableS1_B_antibody_SNPs.csv"), row.names = FALSE)
for (nm in names(IDS)) { if (!is.null(outB[[nm]])) write.csv(outB[[nm]], OUT(paste0("B_out_", nm, ".csv")), row.names = FALSE) }

})

###############################################################################
# C. Reverse MR: EC (larger outcome) -> antibody level (continuous, SD units; beta reported, not OR)
###############################################################################
step("C", {
cat("\n################ C. Reverse MR ################\n")
tryCatch({
  xR <- extract_instruments(ID_LARGE, p1 = 5e-8, clump = TRUE, r2 = 0.001, kb = 10000)
  cat("Number of reverse-MR instruments =", nrow(xR), "(expected: 16)\n")
  dR <- harm(xR, get_out(xR$SNP, ID_AB, TRUE), 2)
  rR <- est(dR, "EC large (12,906)", "Antibody MOMP D (SD)", "C reverse MR: all available SNPs (primary)", primary = TRUE, is_bin = FALSE)
  addres(rR); cat("Reverse MR beta =", round(rR$b, 4), " SE =", round(rR$se, 4), " P =", round(rR$p, 4), " nsnp =", rR$nsnp, "\n")
  write.csv(dR, OUT("C_reverse_harmonised.csv"), row.names = FALSE)
  # Standard Steiger use: drop only SNPs with direction FALSE and P < 0.05 (sensitivity analysis)
  tryCatch({
    dS <- dR; dS$samplesize.outcome <- N_EXPO
    st <- steiger_filtering(dS); print(st[, c("SNP", "rsq.exposure", "rsq.outcome", "steiger_dir", "steiger_pval")])
    drop <- st$SNP[st$steiger_dir %in% FALSE & st$steiger_pval < 0.05]
    cat("SNPs to drop under standard Steiger filtering:", length(drop), "\n")
    addres(est(dR[!(dR$SNP %in% drop), ], "EC large (12,906)", "Antibody MOMP D (SD)",
               "C reverse MR: standard Steiger filtering (sensitivity)", is_bin = FALSE))
  }, error = function(e) cat("Reverse Steiger skipped: ", conditionMessage(e), "\n"))
}, error = function(e) cat("Reverse MR skipped: ", conditionMessage(e), "\n"))

})

###############################################################################
# D. Power (approximation; R2 values are assumed because the true R2 is likely lower than the apparent R2 of the instruments)
###############################################################################
step("D", {
cat("\n################ D. Power ################\n")
getn <- function(o, col) if (is.null(o) || !(col %in% names(o))) NA else suppressWarnings(max(o[[col]], na.rm = TRUE))
pw <- list()
for (nm in c("small", "large", "endo", "nonendo")) {
  o <- if (!is.null(outB[[nm]])) outB[[nm]] else outA[[nm]]
  n <- getn(o, "samplesize.outcome"); nc <- unname(NCASE[nm])      # total sample size from the data, number of cases from NCASE
  if (is.finite(nc) && is.finite(n)) { t <- power_tab(c(0.02, 0.05, 0.10), c(1.05, 1.1, 1.2), nc, n)
    t$outcome <- unname(OUTLAB[nm]); t$ncase <- nc; t$N <- n; pw[[nm]] <- t }
  else cat("Power: skipped", nm, "(sample size missing, N =", n, ")\n")
}
if (length(pw) > 0) {
  PW <- do.call(rbind, pw); PW$power <- round(PW$power, 3); print(PW)
  write.csv(PW, OUT("D_power_table.csv"), row.names = FALSE)
} else cat("Power table is empty: sample size is missing for all outcomes\n")
cat("(Cross-check one cell with the mRnd web tool, e.g. larger outcome, R2 = 0.05, OR = 1.1)\n")

})

###############################################################################
# E. Summary: all estimates + multiple testing + z test of differences + forest plot
###############################################################################
cat("\n################ E. Summary ################\n")
FINAL <- do.call(rbind, RES)
prim <- which(FINAL$primary)
FINAL$p_bonf <- NA; FINAL$p_BH <- NA
FINAL$p_bonf[prim] <- p.adjust(FINAL$p[prim], "bonferroni"); FINAL$p_BH[prim] <- p.adjust(FINAL$p[prim], "BH")
cat("Number of primary tests (multiple-testing family):", length(prim), "\n")
FINAL[, c("b", "se", "OR", "OR_lo", "OR_hi", "I2GX")] <- round(FINAL[, c("b", "se", "OR", "OR_lo", "OR_hi", "I2GX")], 4)
FINAL$p <- signif(FINAL$p, 3); FINAL$p_bonf <- signif(FINAL$p_bonf, 3); FINAL$p_BH <- signif(FINAL$p_BH, 3)
write.csv(FINAL, OUT("FINAL_summary.csv"), row.names = FALSE); print(FINAL[, c("exposure", "outcome", "set", "nsnp", "b", "se", "p")])
ZDF <- do.call(rbind, ZD); write.csv(ZDF, OUT("FINAL_zdiff.csv"), row.names = FALSE); print(ZDF)

tryCatch({
  library(ggplot2)
  fp <- FINAL[FINAL$primary & !is.na(FINAL$OR), ]
  fp$lab <- paste0(fp$outcome, " | ", fp$set, " (n=", fp$nsnp, ")")
  g <- ggplot(fp, aes(x = OR, y = reorder(lab, -seq_len(nrow(fp))))) +
    geom_vline(xintercept = 1, linetype = 2) + geom_point() +
    geom_errorbar(aes(xmin = OR_lo, xmax = OR_hi), orientation = "y", width = 0.2) +
    facet_wrap(~exposure, scales = "free", ncol = 1) +
    labs(x = "OR per unit (reinfection: log-odds; antibody: SD), 95% CI", y = NULL) + theme_bw()
  ggsave(OUT("FINAL_forest.png"), g, width = 9, height = 7, dpi = 300)
}, error = function(e) cat("Forest plot skipped: ", conditionMessage(e), "\n"))

writeLines(capture.output(sessionInfo()), OUT("sessionInfo.txt"))
cat("\nDone. Results are in final_out/ (FINAL_summary.csv, FINAL_zdiff.csv, D_power_table.csv, TableS1_*.csv).\n")
