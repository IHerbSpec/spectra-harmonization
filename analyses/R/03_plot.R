ids <- sort(unique(ctx$metadata$dataset_id))
colors <- setNames(hcl.colors(length(ids),palette="Dark 3"),ids)
write_csv(data.frame(dataset_id=ids,color=unname(colors)),file.path(ctx$out,"tables","dataset_colors.csv"))
fig <- file.path(ctx$out,"figures")
for (type in c("reflectance","normalized")) {
  s <- ctx$summaries[[type]]
  label <- if(type=="reflectance") "Reflectance (fraction)" else "Vector-normalized value"
  for (group in config$material_groups) {
    prefix <- if(group=="pnt") "04" else if(group=="tvk") "05" else "material"
    plot_mean_pages(s[s$material_group==group,,drop=FALSE],
      file.path(fig,paste0(prefix,"_",safe_name(group),"_",type)),colors,label)
  }
  prefix <- if(type=="reflectance") "06" else "10"
  plot_pca(ctx$pca[[type]],file.path(fig,paste0(prefix,"_pca_",type)),colors,
           paste(config$pca_material_group,if(type=="reflectance") "reflectance" else "vector normalized"))
}
for (mode in names(ctx$matched)) {
  ii <- ctx$matched[[mode]]
  prefix <- if(mode=="instruments") "07_matched_instruments" else "08_matched_svc_optics"
  if (!length(ii)) {
    plot_note(file.path(fig,prefix),"No eligible matches. Supply physical_sample_id for confirmed shared targets; optics comparisons require the same instrument_id.")
  } else {
    s <- spectral_summary(ctx$reflectance[ii,,drop=FALSE],ctx$metadata[ii,,drop=FALSE],ctx$wavelengths)
    plot_mean_pages(s,file.path(fig,prefix),colors)
  }
}
plot_differences(ctx$pairs$reflectance$differences,file.path(fig,"paired_mean_differences"))
writeLines(c(
  "# Characterization run", "",
  paste("Included datasets:",length(unique(ctx$metadata$dataset_id))),
  paste("Spectra:",nrow(ctx$metadata)),paste("Output bands:",length(ctx$wavelengths)),"",
  "Figures 04/05: mean and 5th-95th percentile bands, separated by kit/material/physical sample.",
  "Figures 06/10: pnt PCA before and after vector normalization. PCA centers bands without unit-variance scaling.",
  "Figures 07/08: confirmed physical-target matches across instrument families or SVC optics.",
  "Paired differences and metrics compare dataset means on the same physical target.",
  "No constant-offset tests or transfer models are fitted.","",
  "Inspect tables/import_qc.csv and dataset_inventory.csv before interpreting the plots.",
  "Use processed/metadata.csv to recover every sample's source and measurement setup.",
  "*_SKIPPED.txt explains unavailable analyses. Raw files are never modified.",
  "Percentile bands show scan variability, not confidence intervals or independent replication."
),file.path(ctx$out,"README.md"))
writeLines("COMPLETE",file.path(ctx$out,"RUN_STATUS.txt"))
