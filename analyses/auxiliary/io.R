# Shared input contract. One manifest row is one dataset/instrument/optic setup.
resolve_path <- function(path, root) {
  if (grepl("^(/|[A-Za-z]:[/\\\\]|~)",path)) path.expand(path) else file.path(root,path)
}
read_csv <- function(path) read.csv(path, check.names=FALSE, stringsAsFactors=FALSE, colClasses="character", na.strings=c("", "NA"))
write_csv <- function(x,path) write.csv(x,path,row.names=FALSE,na="")
require_columns <- function(x, cols, label) {
  if (!all(cols %in% names(x))) stop(label, " missing: ", paste(setdiff(cols,names(x)),collapse=", "))
  if (anyDuplicated(names(x))) stop(label," has duplicate column names")
}
required_text <- function(x, cols, label) {
  require_columns(x,cols,label)
  for (n in cols) if (any(is.na(x[[n]]) | !nzchar(trimws(x[[n]])))) stop(label,": blank ",n)
}
parse_include <- function(x) {
  if (anyNA(x) || any(!tolower(x) %in% c("true","false"))) stop("include must be TRUE or FALSE")
  tolower(x)=="true"
}
manifest_read <- function(path) {
  x <- read_csv(path)
  required_text(x,c("dataset_id","format","spectra_path","metadata_path","institution",
                    "instrument_family","instrument_model","instrument_id","optic_id",
                    "reflectance_units","include"),"Manifest")
  if (!nrow(x)) stop("Manifest has no datasets")
  if (anyDuplicated(x$dataset_id) || any(!grepl("^[A-Za-z0-9_-]+$",x$dataset_id))) stop("dataset_id must be unique and contain only letters, numbers, hyphens, underscores")
  if (any(!x$format %in% c("wide_csv","long_csv","native"))) stop("Unsupported format")
  if (any(!x$reflectance_units %in% c("fraction","percent"))) stop("Set reflectance_units to fraction or percent explicitly")
  x$include <- parse_include(x$include)
  if (!any(x$include)) stop("No included datasets")
  if (!"exclusion_reason" %in% names(x)) x$exclusion_reason <- ""
  if (any(!x$include & (is.na(x$exclusion_reason)|!nzchar(x$exclusion_reason)))) stop("Excluded datasets need exclusion_reason")
  x
}

metadata_read <- function(path) {
  m <- read_csv(path)
  # Aliases accept metadata exported by this repository's Shiny QC app.
  aliases <- c(spectrum_id="fullFilename",material_id="targetClass",kit_id="kitNumber")
  for (n in names(aliases)) if (!n %in% names(m) && aliases[[n]] %in% names(m)) m[[n]] <- m[[aliases[[n]]]]
  required_text(m,c("spectrum_id","material_id","kit_id"),"Metadata")
  if (!nrow(m)) stop("Metadata has no spectra")
  if (anyDuplicated(m$spectrum_id)) stop("Duplicate spectrum_id in metadata")
  m$material_id <- sub("^TC","",m$material_id)
  if (!"material_group" %in% names(m)) m$material_group <- ifelse(grepl("^pnt",m$material_id),"pnt",m$material_id)
  required_text(m,"material_group","Metadata")
  if (!"physical_sample_id" %in% names(m)) m$physical_sample_id <- NA_character_
  # Missing physical IDs remain unmatched, rather than implying identical targets.
  m$physical_sample_id[!is.na(m$physical_sample_id) & !nzchar(trimws(m$physical_sample_id))] <- NA_character_
  reserved <- c("dataset_id","sample_key","institution","instrument_family","instrument_model","instrument_id","optic_id","reflectance_units")
  if (any(reserved %in% names(m))) stop("Metadata contains manifest-reserved fields: ",paste(intersect(reserved,names(m)),collapse=", "))
  m
}

read_wide <- function(path) {
  d <- read_csv(path)
  required_text(d,"spectrum_id","Wide spectra")
  cols <- grepl("^(X)?[0-9]+(\\.[0-9]+)?$",names(d))
  if (sum(cols)<2L) stop("Wide CSV needs numeric wavelength headers, e.g. 400 or X400")
  if (sum(!cols)!=1L) stop("Wide spectra CSV must contain only spectrum_id and wavelength columns; keep metadata separately")
  w <- as.numeric(sub("^X","",names(d)[cols]))
  x <- as.matrix(data.frame(lapply(d[cols],as.numeric),check.names=FALSE))
  ord <- order(w)
  list(ids=d$spectrum_id,wavelengths=w[ord],spectra=x[,ord,drop=FALSE],files=path)
}

read_long <- function(path) {
  d <- read_csv(path)
  required_text(d,c("spectrum_id","wavelength_nm","reflectance"),"Long spectra")
  d$wavelength_nm <- as.numeric(d$wavelength_nm)
  d$reflectance <- as.numeric(d$reflectance)
  if (any(!is.finite(d$wavelength_nm)) || any(!is.finite(d$reflectance))) stop("Non-numeric long spectra")
  if (anyDuplicated(d[c("spectrum_id","wavelength_nm")])) stop("Duplicate sample/wavelength in long CSV")
  ids <- unique(d$spectrum_id); w <- sort(unique(d$wavelength_nm))
  x <- matrix(NA_real_,length(ids),length(w))
  x[cbind(match(d$spectrum_id,ids),match(d$wavelength_nm,w))] <- d$reflectance
  if (anyNA(x)) stop("Long CSV has inconsistent wavelength grids; split datasets by grid")
  list(ids=ids,wavelengths=w,spectra=x,files=path)
}

read_native <- function(path, metadata) {
  if (!requireNamespace("spectrolab",quietly=TRUE)) stop("Install spectrolab to import native files: install.packages('spectrolab')")
  # Read one explicitly named file at a time; never infer metadata order from a
  # directory listing or from sample labels returned by a native-file reader.
  files <- file.path(path,metadata$spectrum_id)
  if (any(basename(metadata$spectrum_id)!=metadata$spectrum_id)) stop("Native spectrum_id must be a filename, not a path")
  if (any(!file.exists(files))) stop("Missing native files: ",paste(basename(files[!file.exists(files)]),collapse=", "))
  values <- lapply(files,function(f) {
    ext <- tolower(tools::file_ext(f))
    if (!ext %in% c("sed","sig","asd","txt")) stop("Unsupported native extension: ",ext)
    obj <- spectrolab::read_spectra(path=f,format=ext,type="target_reflectance",extract_metadata=FALSE)
    if (!inherits(obj,"spectra")) stop("Reader did not return a spectra object: ",f)
    mat <- as.matrix(obj)
    if (nrow(mat)!=1L) stop("Expected one spectrum per native file: ",f)
    list(w=as.numeric(spectrolab::bands(obj)),x=mat)
  })
  w <- values[[1]]$w
  if (!all(vapply(values,function(v) identical(v$w,w),logical(1)))) stop("Different grids within native dataset; split manifest rows by grid")
  ord <- order(w)
  list(ids=metadata$spectrum_id,wavelengths=w[ord],spectra=do.call(rbind,lapply(values,`[[`,"x"))[,ord,drop=FALSE],files=files)
}

load_dataset <- function(row,root,config) {
  mpath <- resolve_path(row$metadata_path,root)
  spath <- resolve_path(row$spectra_path,root)
  m <- metadata_read(mpath)
  d <- switch(row$format,wide_csv=read_wide(spath),long_csv=read_long(spath),native=read_native(spath,m))
  if (anyDuplicated(d$ids) || !setequal(d$ids,m$spectrum_id)) stop(row$dataset_id,": spectra and metadata IDs must match exactly and uniquely")
  m <- m[match(d$ids,m$spectrum_id),,drop=FALSE]
  for (n in c("dataset_id","institution","instrument_family","instrument_model","instrument_id","optic_id","reflectance_units")) m[[n]] <- row[[n]]
  m$sample_key <- paste(row$dataset_id,m$spectrum_id,sep="::")
  x <- d$spectra
  if (row$reflectance_units=="percent") x <- x/100
  if (any(!is.finite(x))) stop(row$dataset_id,": missing/non-numeric reflectance; resolve before processing")
  # Preserve out-of-range values for characterization, and report them.
  qc <- data.frame(sample_key=m$sample_key,source_min_nm=min(d$wavelengths),source_max_nm=max(d$wavelengths),
                   source_max_gap_nm=max(diff(d$wavelengths)),source_n_bands=length(d$wavelengths),
                   native_min_reflectance=apply(x,1,min),native_max_reflectance=apply(x,1,max),
                   n_outside_0_1=rowSums(x<0 | x>1))
  target <- seq(config$wavelength_min_nm,config$wavelength_max_nm,by=config$spacing_nm)
  y <- gaussian_resample(x,d$wavelengths,target,config$target_fwhm_nm,config$max_source_gap_nm)
  list(metadata=m,spectra=y,qc=qc,files=c(mpath,d$files))
}

bind_metadata <- function(items) {
  cols <- unique(unlist(lapply(items,names)))
  do.call(rbind,lapply(items,function(x) {x[setdiff(cols,names(x))] <- NA_character_;x[,cols,drop=FALSE]}))
}
