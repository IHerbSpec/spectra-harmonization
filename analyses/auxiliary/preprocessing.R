# Gaussian target-band resampling. See docs/methods.md for the distinction
# between target bandpass, output spacing, and native instrument resolution.
validate_config <- function(x) {
  required <- c("manifest", "output_dir", "wavelength_min_nm", "wavelength_max_nm",
                "spacing_nm", "target_fwhm_nm", "max_source_gap_nm",
                "min_replicates", "material_groups", "pca_material_group", "overwrite")
  if (!all(required %in% names(x))) stop("Missing config fields: ", paste(setdiff(required, names(x)), collapse=", "))
  for (n in c("wavelength_min_nm", "wavelength_max_nm", "spacing_nm", "target_fwhm_nm",
              "max_source_gap_nm", "min_replicates"))
    if (length(x[[n]]) != 1L || !is.finite(x[[n]]) || x[[n]] <= 0) stop("Invalid config: ", n)
  if (x$wavelength_max_nm <= x$wavelength_min_nm) stop("Invalid wavelength range")
  intervals <- (x$wavelength_max_nm-x$wavelength_min_nm)/x$spacing_nm
  if (abs(intervals-round(intervals)) > 1e-8) stop("Wavelength endpoints must align with spacing")
  if (!is.logical(x$overwrite) || length(x$overwrite)!=1L || is.na(x$overwrite)) stop("overwrite must be TRUE or FALSE")
  x
}

gaussian_resample <- function(spectra, wavelengths, new_wavelengths,
                              fwhm_nm = 5, max_gap_nm = 5) {
  x <- as.matrix(spectra)
  w <- as.numeric(wavelengths)
  target <- as.numeric(new_wavelengths)
  if (!is.numeric(x) || ncol(x)!=length(w) || !nrow(x) || length(w)<2L ||
      any(!is.finite(w)) || any(!is.finite(x))) stop("Non-finite or invalid spectral input")
  if (any(diff(w)<=0)) stop("Source wavelengths must be unique and increasing")
  if (!length(target) || any(!is.finite(target)) || any(diff(target)<=0)) stop("Invalid target wavelengths")
  if (!(length(fwhm_nm) %in% c(1L,length(target))) || any(!is.finite(fwhm_nm)) || any(fwhm_nm<=0)) stop("Invalid FWHM")
  fwhm_nm <- rep(fwhm_nm, length.out=length(target))
  sigma <- fwhm_nm/(2*sqrt(2*log(2)))
  if (any(target-3*sigma < min(w)-1e-8 | target+3*sigma > max(w)+1e-8))
    stop("Insufficient wavelength coverage for complete Gaussian windows. Retain native bands beyond the final limits (about 394-2406 nm for FWHM=5).")
  # Check only the working range, including kernel margins.
  use <- w >= max(min(w),min(target-3*sigma)-max_gap_nm) &
    w <= min(max(w),max(target+3*sigma)+max_gap_nm)
  if (any(diff(w[use])>max_gap_nm+1e-8)) stop("Source wavelength gap exceeds max_source_gap_nm; do not bridge missing detector regions")
  # Discrete channel weights reproduce CWT::resampling_FWHM's weighting.
  # Retain all native bands for numerical agreement and endpoint support.
  weights <- vapply(seq_along(target), function(j) {
    z <- exp(-0.5*((w-target[j])/sigma[j])^2)
    z/sum(z)
  }, numeric(length(w)))
  out <- x %*% weights
  colnames(out) <- as.character(target)
  out
}

vector_normalize <- function(x) {
  norms <- sqrt(rowSums(x*x))
  if (any(!is.finite(norms) | norms<=0)) stop("Cannot vector-normalize zero or non-finite spectrum")
  list(spectra=x/norms, norms=norms)
}
