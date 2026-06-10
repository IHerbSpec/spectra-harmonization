# ============================================================
# qc_spectra_files.R
#
# Quality control for spectral files generated from
# IHerbSpec harmonization kits.
#
# Recommended use:
#   Open spectra-harmonization.Rproj, then run this script.
#   All paths are relative to the repository root.
#
# Workflow:
#   1. Edit User settings below.
#   2. Run Steps 1–3 to generate QC plots.
#   3. Review plots and edit files_to_delete.csv as needed.
# ============================================================

library(spectrolab)
library(tidyverse)
library(ggrepel)

source(file.path("scripts", "qc_functions.R"))

# ============================================================
# User settings  ← edit these
# ============================================================

min_expected_files <- 5        # Minimum measurements expected per material

# ============================================================
# Paths  (relative to repository root — do not edit)
# ============================================================

spectra_dir <- "raw_data_files"   # folder containing .sig/.sed/.asd/.txt files
out_dir     <- "outputs"
qc_plot_dir <- file.path(out_dir, "qc_plots")

# ============================================================
# Check required directories and files
# ============================================================

required_dirs <- c(spectra_dir, out_dir, qc_plot_dir)
missing_dirs  <- required_dirs[!dir.exists(required_dirs)]

if (length(missing_dirs) > 0) {
  stop(
    "Missing required directories:\n",
    paste(missing_dirs, collapse = "\n"),
    "\nThese should be present in the cloned repository."
  )
}


# ============================================================
# Step 1: Parse filenames
# ============================================================

file_inventory <- list_spectral_files(spectra_dir)
cat("Spectral files found:", nrow(file_inventory), "\n")

parsed        <- parse_filenames(file_inventory)
valid_files   <- parsed$valid
bad_filenames <- parsed$bad

write_csv(bad_filenames, file.path(out_dir, "bad_filenames.csv"))

if (nrow(bad_filenames) > 0) {
  warning("Some files do not match expected filename conventions. See bad_filenames.csv")
  print(bad_filenames, n = Inf)
}

if (nrow(valid_files) == 0) stop("No files matched the expected filename conventions.")

# ============================================================
# Step 2: Count files per material
# ============================================================

file_counts       <- check_file_counts(valid_files, min_expected_files)
flagged_materials <- file_counts |> filter(n_files < min_expected_files)

write_csv(file_counts, file.path(out_dir, "file_counts_by_material.csv"))

cat("\nFile counts by material:\n")
print(file_counts, n = Inf)

cat("\nFlagged materials (fewer than", min_expected_files, "files):\n")
print(flagged_materials, n = Inf)

missing_materials <- check_missing_materials(valid_files)
if (nrow(missing_materials) > 0) {
  write_csv(missing_materials, file.path(out_dir, "missing_materials_warning.csv"))
  warning(
    nrow(missing_materials), " of 19 expected materials not found: ",
    paste(missing_materials$missing_material, collapse = ", "),
    ". See missing_materials_warning.csv"
  )
}

# ============================================================
# Step 3: Read spectral files and generate QC plots
# ============================================================

spec         <- read_spectral_files(valid_files)
spectra_long <- spectra_to_long(spec, valid_files)

pdf_path <- save_qc_plots(spectra_long, qc_plot_dir, out_dir)

cat("\nQC plots saved.\n")
cat("Multipage PDF:   ", pdf_path, "\n")
cat("Individual PNGs: ", qc_plot_dir, "\n")

# ============================================================
# Done
# ============================================================

cat("\nQC complete. Key outputs:\n")
cat(" ", file.path(out_dir, "file_counts_by_material.csv"), "\n")
cat(" ", file.path(out_dir, "bad_filenames.csv"), "\n")
cat(" ", pdf_path, "\n")
