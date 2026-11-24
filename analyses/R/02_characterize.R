# Summaries retain kit and material, so different kits are not pooled silently.
ctx$summaries <- list();ctx$pca <- list();ctx$pairs <- list()
for (type in c("reflectance","normalized")) {
  x <- ctx[[type]]
  ctx$summaries[[type]] <- spectral_summary(x,ctx$metadata,ctx$wavelengths)
  write_csv(ctx$summaries[[type]],file.path(ctx$out,"tables",paste0("summary_",type,".csv")))
  ii <- which(ctx$metadata$material_group==config$pca_material_group)
  fit <- fit_pca(x[ii,,drop=FALSE],ctx$metadata[ii,,drop=FALSE])
  ctx$pca[[type]] <- fit
  if (!is.null(fit$skipped)) {
    writeLines(fit$skipped,file.path(ctx$out,"tables",paste0("pca_",type,"_SKIPPED.txt")))
  } else {
    for (part in c("scores","loadings","variance")) write_csv(fit[[part]],file.path(ctx$out,"tables",paste0("pca_",type,"_",part,".csv")))
    saveRDS(fit$model,file.path(ctx$out,"processed",paste0("pca_",type,".rds")))
  }
  ctx$pairs[[type]] <- pairwise_means(x,ctx$metadata,ctx$wavelengths)
  for (part in c("metrics","differences")) {
    d <- ctx$pairs[[type]][[part]]
    if (nrow(d)) write_csv(d,file.path(ctx$out,"tables",paste0("paired_",type,"_",part,".csv")))
  }
}
inventory_cols <- c("dataset_id","institution","instrument_family","instrument_model","instrument_id","optic_id","kit_id","material_id","material_group")
ctx$inventory <- do.call(rbind,lapply(group_indices(ctx$metadata,inventory_cols),function(ii) {
  cbind(ctx$metadata[ii[1],inventory_cols,drop=FALSE],n_spectra=length(ii),below_min_replicates=length(ii)<config$min_replicates)
}))
write_csv(ctx$inventory,file.path(ctx$out,"tables","dataset_inventory.csv"))
write_csv(data.frame(sample_key=ctx$metadata$sample_key,mean_reflectance=rowMeans(ctx$reflectance),
                     vector_norm=sqrt(rowSums(ctx$reflectance^2))),file.path(ctx$out,"tables","magnitude_per_spectrum.csv"))
ctx$matched <- lapply(c("instruments","svc_optics"),function(mode) matched_rows(ctx$metadata,mode))
names(ctx$matched) <- c("instruments","svc_optics")
for (mode in names(ctx$matched)) {
  ii <- ctx$matched[[mode]]
  write_csv(ctx$metadata[ii,,drop=FALSE],file.path(ctx$out,"tables",paste0("matched_",mode,"_samples.csv")))
}
