# All paths are relative to analyses/, unless absolute.
# Copy this file to config/local.R for your own datasets.
config <- list(
  manifest = "config/datasets.example.csv",
  output_dir = "outputs/example",
  wavelength_min_nm = 400,
  wavelength_max_nm = 2400,
  spacing_nm = 5,
  target_fwhm_nm = 5,   # Gaussian kernel width, independent of spacing
  max_source_gap_nm = 5,
  min_replicates = 5,
  material_groups = c("pnt", "tvk"),
  pca_material_group = "pnt",
  overwrite = FALSE
)
