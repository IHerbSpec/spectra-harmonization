# Import, validate, Gaussian-resample, then retain the requested band range.
ctx$manifest_path <- resolve_path(config$manifest,root)
ctx$manifest <- manifest_read(ctx$manifest_path)
ctx$out <- resolve_path(config$output_dir,root)
if (dir.exists(ctx$out) && length(list.files(ctx$out,all.files=TRUE,no..=TRUE))) {
  if (!config$overwrite) stop("Output directory is not empty. Select a new output_dir or explicitly set overwrite=TRUE.")
  # Only remove generated subfolders inside the declared output directory.
  for (d in c("processed","tables","figures","provenance")) unlink(file.path(ctx$out,d),recursive=TRUE)
}
for (d in c("processed","tables","figures","provenance")) dir.create(file.path(ctx$out,d),recursive=TRUE,showWarnings=FALSE)
writeLines("RUNNING",file.path(ctx$out,"RUN_STATUS.txt"))
write_csv(ctx$manifest,file.path(ctx$out,"provenance","dataset_manifest.csv"))
active <- ctx$manifest[ctx$manifest$include,,drop=FALSE]
datasets <- lapply(seq_len(nrow(active)),function(i) {
  message("  Importing ",active$dataset_id[i])
  tryCatch(load_dataset(active[i,,drop=FALSE],root,config),error=function(e) stop(active$dataset_id[i],": ",conditionMessage(e),call.=FALSE))
})
ctx$metadata <- bind_metadata(lapply(datasets,`[[`,"metadata"))
ctx$reflectance <- do.call(rbind,lapply(datasets,`[[`,"spectra"))
ctx$wavelengths <- as.numeric(colnames(ctx$reflectance))
if (anyDuplicated(ctx$metadata$sample_key)) stop("Sample keys are not unique")
# The same physical ID must not identify distinct kits, materials, or surfaces.
ids <- ctx$metadata[!is.na(ctx$metadata$physical_sample_id),c("physical_sample_id","kit_id","material_id")]
if (nrow(ids) && anyDuplicated(unique(ids)$physical_sample_id)) stop("A physical_sample_id maps to multiple kits or materials")
instruments <- unique(ctx$metadata[c("instrument_id","instrument_family","instrument_model","institution")])
if (anyDuplicated(instruments$instrument_id)) stop("instrument_id must identify a single physical instrument and institution")
vn <- vector_normalize(ctx$reflectance)
ctx$normalized <- vn$spectra
ctx$qc <- do.call(rbind,lapply(datasets,`[[`,"qc"))
ctx$qc$vector_norm <- vn$norms
if (any(ctx$qc$n_outside_0_1>0)) warning("Reflectance outside [0,1] retained. Inspect tables/import_qc.csv; units may need checking.")
write_csv(ctx$qc,file.path(ctx$out,"tables","import_qc.csv"))
write_csv(ctx$metadata,file.path(ctx$out,"processed","metadata.csv"))
for (type in c("reflectance","normalized")) {
  mat <- ctx[[type]]
  write_csv(data.frame(sample_key=ctx$metadata$sample_key,mat,check.names=FALSE),file.path(ctx$out,"processed",paste0(type,"_5nm.csv")))
}
saveRDS(list(metadata=ctx$metadata,wavelengths=ctx$wavelengths,reflectance=ctx$reflectance,
             normalized=ctx$normalized,config=config),file.path(ctx$out,"processed","spectra.rds"))
inputs <- unique(c(ctx$manifest_path,config_file,unlist(lapply(datasets,`[[`,"files"))))
write_csv(data.frame(path=inputs,md5=unname(tools::md5sum(inputs))),file.path(ctx$out,"provenance","input_checksums.csv"))
scripts <- c(list.files(file.path(root,"R"),full.names=TRUE,pattern="[.]R$"),
  list.files(file.path(root,"auxiliary"),full.names=TRUE,pattern="[.]R$"),
  list.files(file.path(root,"plotting"),full.names=TRUE,pattern="[.]R$"),file.path(root,"run_all.R"))
write_csv(data.frame(path=scripts,md5=unname(tools::md5sum(scripts))),file.path(ctx$out,"provenance","code_checksums.csv"))
dput(config,file=file.path(ctx$out,"provenance","config.R"))
writeLines(capture.output(sessionInfo()),file.path(ctx$out,"provenance","sessionInfo.txt"))
writeLines(c(paste("UTC:",format(Sys.time(),tz="UTC",usetz=TRUE)),
             "Method: discrete Gaussian target-band integration; no native-resolution matching",
             paste("Target grid:",paste(range(ctx$wavelengths),collapse=" to "),"nm"),
             paste("Spacing:",config$spacing_nm,"nm; target Gaussian FWHM:",config$target_fwhm_nm,"nm")),
           file.path(ctx$out,"provenance","run.txt"))
