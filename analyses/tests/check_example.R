local({
  out <- "analyses/outputs/example"
  stopifnot(identical(readLines(file.path(out,"RUN_STATUS.txt")),"COMPLETE"))
  d <- readRDS(file.path(out,"processed/spectra.rds"))
  stopifnot(nrow(d$metadata)==120L,ncol(d$reflectance)==401L,
            identical(d$wavelengths,seq(400,2400,5)),
            max(abs(sqrt(rowSums(d$normalized^2))-1))<1e-12)
  expected <- c("figures/04_pnt_reflectance.pdf","figures/05_tvk_reflectance.pdf",
    "figures/06_pca_reflectance.png","figures/10_pca_normalized.png",
    "figures/07_matched_instruments.pdf","figures/08_matched_svc_optics.pdf",
    "tables/paired_reflectance_metrics.csv","tables/summary_normalized.csv",
    "provenance/input_checksums.csv","provenance/code_checksums.csv")
  stopifnot(all(file.exists(file.path(out,expected))))
  paired <- read.csv(file.path(out,"tables/paired_reflectance_metrics.csv"))
  stopifnot(nrow(paired)==36L,all(is.finite(paired$rmse)))
  message("PASS: complete synthetic workflow, 120 spectra, 401 bands, 36 matched comparisons")
})
