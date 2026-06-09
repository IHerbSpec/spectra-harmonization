# ============================================================
# qc_spectra_files.R
#
# Quality control for spectral files generated from
# IHerbSpec harmonization kits.
#
# This script:
#   1. Checks filename conventions (outputs parsed filename info and problem files)
#   2. Counts files per material
#   3. Flags materials with fewer than the expected number of files (outputs list to csv)
#   4. Reads spectral files with spectrolab (.sed, .sig, .asd, or .txt)
#   5. Generates labeled QC plots for visual inspection for outlier spectra
#   6. List/flag visually abnormal files in csv
#   7. Exports good files with converted full filename convention
# ============================================================

# ---- packages ----
# install.packages(c("spectrolab", "tidyverse", "stringr", "ggrepel"))

library(spectrolab)
library(tidyverse)
library(stringr)
library(ggrepel)

# ============================================================
# Working directory
# ============================================================
#
# Recommended use:
#   1. Open spectra-harmonization.Rproj in RStudio.
#   2. Run this script from the repository root.
#
# Check your current working directory with:
#
#   getwd()
#
# The script assumes paths are relative to the repository root.
#
# ============================================================
# User settings
# ============================================================

### Revise these:

# Herbarium code for the participating lab (e.g. "HCMO")
herbariumCode <- "HCHUH"

# kitNumber (this will be parsed from full filename if it exists)
kitNumber <- "1"

# Minimum expected number of measurements per targetID + kitNumber.
min_expected_files <- 5

### Run these:

# Folder containing raw spectral files.
spectra_dir <- file.path("data-collection-QC", "raw_data_files")

# Output directory.
out_dir <- file.path("data-collection-QC", "outputs")

# QC plot directory.
qc_plot_dir <- file.path(out_dir, "qc_plots")

# Folder where good files will be copied using full harmonization filenames.
good_files_dir <- file.path(out_dir, "good_files_full_filenames")

# File where users list bad files to exclude after visual QC.
delete_list_csv <- file.path(out_dir, "files_to_delete.csv")

# ============================================================
# Check for required repo folders/files
# ============================================================
#
# These folders/files should already exist in the cloned repository.
# The script does not create them.

required_dirs <- c(
  spectra_dir,
  out_dir,
  qc_plot_dir,
  good_files_dir
)

missing_dirs <- required_dirs[!dir.exists(required_dirs)]

if (length(missing_dirs) > 0) {
  stop(
    "The following required directories are missing:\n",
    paste(missing_dirs, collapse = "\n"),
    "\n\nThese should be present in the repository before running the script."
  )
}

if (!file.exists(delete_list_csv)) {
  stop(
    "Could not find files_to_delete.csv at:\n",
    delete_list_csv,
    "\n\nThis file should already exist and contain columns: filename_with_extension, reason"
  )
}

# ============================================================
# File formats and filename parsing
# ============================================================
#
# Accepted file extensions:
#   .sed, .sig, .asd, .txt
#
# Supported filename structures:
#
# 1. Simple local QC filenames:
#      targetID_IDX.ext
#
#    Example:
#      fab2_0002.sig
#      pap1_0010.sed
#      tyvek_0005.sig
#
#    In this case, kitNumber is not included in the filename.
#    The script assigns the kitNumber entered in the User settings.
#
# 2. Full harmonization filenames for sharing:
#      PIdataHarmonization2026_HCHUH_targetID_KitNumber_IDX.ext
#
#    Example:
#      PIdataHarmonization2026_HCHUH_fab2_1_0002.sig
#      PIdataHarmonization2026_HCMO_tyvek_1_0005.sed
#
# The projectID, PIdataHarmonization2026, is fixed for all files
# from all participants.
#
# The herbariumCode, HC..., is fixed per institution.
#
# The script uses targetID + kitNumber as the material_id for QC.
# IDX identifies the replicate measurement file. IDX can be any
# numeric value and does not need to be sequential or 01-05.
#
# QC checks whether each targetID + kitNumber has at least
# min_expected_files files, regardless of IDX values.
#
# ============================================================

valid_extensions <- c("sed", "sig", "asd", "txt")

spectral_files <- list.files(
  spectra_dir,
  pattern = paste0("\\.(", paste(valid_extensions, collapse = "|"), ")$"),
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(spectral_files) == 0) {
  stop(
    "No spectral files found in: ", spectra_dir,
    "\nExpected extensions: .sed, .sig, .asd, .txt"
  )
}

file_inventory <- tibble(
  path = spectral_files,
  file = basename(spectral_files),
  extension = str_to_lower(str_extract(file, "[^.]+$")),
  filename_no_ext = str_remove(file, "\\.[^.]+$")
)

cat("Number of spectral files found:", nrow(file_inventory), "\n")

# ============================================================
# Parse filenames
# ============================================================
#
# Full harmonization filename:
#   PIdataHarmonization2026_HCHUH_targetID_KitNumber_IDX
#
# Simple local filename:
#   targetID_IDX
#
# targetID is allowed to contain underscores.
# The final numeric field is always IDX.
# In full harmonization names, the numeric field before IDX is kitNumber.
#
# ============================================================

shared_regex <- "^(PIdataHarmonization2026)_([A-Za-z0-9]+)_(.+)_(\\d+)_(\\d+)$"
local_regex  <- "^(.+)_(\\d+)$"

file_inventory <- file_inventory %>%
  mutate(
    is_shared_name = str_detect(filename_no_ext, shared_regex),
    is_local_name = str_detect(filename_no_ext, local_regex) & !is_shared_name,
    filename_ok = is_shared_name | is_local_name
  ) %>%
  mutate(
    projectID_parsed = if_else(
      is_shared_name,
      str_match(filename_no_ext, shared_regex)[, 2],
      NA_character_
    ),
    
    herbariumCode_parsed = if_else(
      is_shared_name,
      str_match(filename_no_ext, shared_regex)[, 3],
      NA_character_
    ),
    
    targetID_shared = if_else(
      is_shared_name,
      str_match(filename_no_ext, shared_regex)[, 4],
      NA_character_
    ),
    
    kitNumber_shared = if_else(
      is_shared_name,
      str_match(filename_no_ext, shared_regex)[, 5],
      NA_character_
    ),
    
    IDX_shared = if_else(
      is_shared_name,
      str_match(filename_no_ext, shared_regex)[, 6],
      NA_character_
    ),
    
    targetID_local = if_else(
      is_local_name,
      str_match(filename_no_ext, local_regex)[, 2],
      NA_character_
    ),
    
    kitNumber_local = if_else(
      is_local_name,
      kitNumber,
      NA_character_
    ),
    
    IDX_local = if_else(
      is_local_name,
      str_match(filename_no_ext, local_regex)[, 3],
      NA_character_
    )
  ) %>%
  mutate(
    targetID = coalesce(targetID_shared, targetID_local),
    kitNumber_parsed_or_entered = coalesce(kitNumber_shared, kitNumber_local),
    IDX = coalesce(IDX_shared, IDX_local),
    
    kitNumber_int = as.integer(kitNumber_parsed_or_entered),
    IDX_int = as.integer(IDX),
    
    material_id = paste(targetID, kitNumber_parsed_or_entered, sep = "_"),
    
    filename_type = case_when(
      is_shared_name ~ "full_harmonization_filename",
      is_local_name ~ "simple_local_filename",
      TRUE ~ "unrecognized_filename"
    )
  )

# ============================================================
# Identify bad filenames
# ============================================================

bad_filenames <- file_inventory %>%
  filter(
    !filename_ok |
      is.na(targetID) |
      is.na(kitNumber_int) |
      is.na(IDX_int)
  )

write_csv(
  bad_filenames,
  file.path(out_dir, "bad_filenames.csv")
)

if (nrow(bad_filenames) > 0) {
  warning("Some files do not match the expected filename conventions. See bad_filenames.csv")
  print(bad_filenames, n = Inf)
}

valid_files <- file_inventory %>%
  filter(
    filename_ok,
    !is.na(targetID),
    !is.na(kitNumber_int),
    !is.na(IDX_int)
  )

if (nrow(valid_files) == 0) {
  stop("No files matched the expected filename conventions.")
}

write_csv(
  valid_files,
  file.path(out_dir, "parsed_file_inventory.csv")
)

cat("\nParsed filename summary:\n")
valid_files %>%
  count(filename_type, extension, name = "n_files") %>%
  print(n = Inf)

# Optional warning if full filenames use a different kitNumber than the user setting.
# This does not stop the script; full filenames retain their parsed kitNumber.

shared_kit_mismatch <- valid_files %>%
  filter(
    filename_type == "full_harmonization_filename",
    kitNumber_parsed_or_entered != kitNumber
  )

if (nrow(shared_kit_mismatch) > 0) {
  warning(
    "Some full harmonization filenames have a kitNumber different from the user-entered kitNumber. ",
    "Parsed kitNumber from the filename will be used for those files."
  )
  
  shared_kit_mismatch %>%
    distinct(file, targetID, kitNumber_parsed_or_entered, IDX) %>%
    print(n = Inf)
}

# ============================================================
# Count files per material
# ============================================================

file_counts <- valid_files %>%
  count(
    projectID_parsed,
    herbariumCode_parsed,
    targetID,
    kitNumber_parsed_or_entered,
    material_id,
    name = "n_files"
  ) %>%
  mutate(
    qc_flag = if_else(
      n_files < min_expected_files,
      paste0("FLAG: fewer than ", min_expected_files, " files"),
      "OK"
    )
  ) %>%
  arrange(qc_flag, targetID, kitNumber_parsed_or_entered)

flagged_materials <- file_counts %>%
  filter(n_files < min_expected_files)

write_csv(
  file_counts,
  file.path(out_dir, "file_counts_by_material.csv")
)

cat("\nFile counts by material:\n")
print(file_counts, n = Inf)

cat("\nFlagged materials:\n")
print(flagged_materials, n = Inf)


# ============================================================
# Read spectral files by extension
# ============================================================

read_one_extension <- function(files_df, ext) {
  these_files <- files_df %>%
    filter(extension == ext)
  
  if (nrow(these_files) == 0) {
    return(NULL)
  }
  
  cat("Reading", nrow(these_files), paste0(".", ext), "files\n")
  
  read_spectra(
    path = these_files$path,
    format = ext,
    type = "target_reflectance",
    extract_metadata = FALSE
  )
}

spec_list <- map(
  valid_extensions,
  ~ read_one_extension(valid_files, .x)
)

names(spec_list) <- valid_extensions
spec_list <- compact(spec_list)

if (length(spec_list) == 0) {
  stop("No readable spectral files found.")
}

if (length(spec_list) == 1) {
  spec <- spec_list[[1]]
} else {
  spec <- reduce(spec_list, spectrolab::combine)
}

if (is.list(spec) && !inherits(spec, "spectra")) {
  stop(
    "read_spectra() returned a list, likely because files have incompatible wavelength grids. ",
    "Inspect or resample files before plotting."
  )
}

# ============================================================
# Convert spectra to long format
# ============================================================

spec_df <- as.data.frame(
  spec,
  fix_names = "none",
  metadata = FALSE
)

names(spec_df)[1] <- "sample_name"

spectra_long <- spec_df %>%
  pivot_longer(
    cols = -sample_name,
    names_to = "wavelength",
    values_to = "reflectance"
  ) %>%
  mutate(
    wavelength = as.numeric(wavelength),
    file = basename(sample_name)
  ) %>%
  left_join(
    valid_files %>%
      select(
        file,
        extension,
        filename_type,
        projectID_parsed,
        herbariumCode_parsed,
        targetID,
        kitNumber_parsed_or_entered,
        kitNumber_int,
        IDX,
        IDX_int,
        material_id
      ),
    by = "file"
  ) %>%
  filter(!is.na(material_id))

write_csv(
  spectra_long,
  file.path(out_dir, "spectra_long.csv")
)

# ============================================================
# Plot each material with filename labels
# ============================================================

plot_material_qc <- function(material_name) {
  dat <- spectra_long %>%
    filter(material_id == material_name) %>%
    mutate(
      file_label = str_remove(file, paste0("\\.", extension, "$"))
    )
  
  if (nrow(dat) == 0) {
    return(NULL)
  }
  
  label_df <- dat %>%
    group_by(file, file_label) %>%
    filter(wavelength == max(wavelength, na.rm = TRUE)) %>%
    ungroup()
  
  ggplot(
    dat,
    aes(
      x = wavelength,
      y = reflectance,
      group = file,
      color = file
    )
  ) +
    geom_line(linewidth = 0.55, alpha = 0.85, show.legend = FALSE) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = file_label),
      hjust = 0,
      nudge_x = 30,
      direction = "y",
      size = 3.2,
      segment.size = 0.25,
      segment.alpha = 0.5,
      show.legend = FALSE
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.25))) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = paste("QC plot:", material_name),
      subtitle = "Each line is one spectral file; labels identify filenames for visual QC",
      x = "Wavelength (nm)",
      y = "Reflectance"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank()
    )
}

materials <- spectra_long %>%
  distinct(material_id) %>%
  arrange(material_id) %>%
  pull(material_id)

for (mat in materials) {
  p <- plot_material_qc(mat)
  
  safe_mat <- str_replace_all(mat, "[^A-Za-z0-9]+", "_")
  
  ggsave(
    filename = file.path(qc_plot_dir, paste0("QC_", safe_mat, ".png")),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
  )
}

# ============================================================
# Generate multipage PDF with all QC plots
# ============================================================

all_qc_pdf <- file.path(out_dir, "ALL_materials_QC_labeled_spectra.pdf")

pdf(all_qc_pdf, width = 10, height = 6)

for (mat in materials) {
  print(plot_material_qc(mat))
}

dev.off()

# ============================================================
# Review QC plots and list files for exclusion in csv
# ============================================================
#
# At this point, QC plots have been generated.
#
# Users should inspect:
#   data-collection-QC/outputs/ALL_materials_QC_labeled_spectra.pdf
#
# and/or the individual plots in:
#   data-collection-QC/outputs/qc_plots/
#
#
# !! If any spectra look abnormal, ADD  their filenames to: !!
#
#   data-collection-QC/outputs/files_to_delete.csv
#
#
# The CSV must contain only:
#   filename_with_extension, reason (optional)
#
# Example:
#   pnt1_0005.sed,bad NIR
#
# ============================================================

cat("\nQC plots generated.\n")
cat("Review the multipage QC PDF:\n")
cat(all_qc_pdf, "\n")

cat("\nReview individual QC plots in:\n")
cat(qc_plot_dir, "\n")

cat("\nIf any spectra look abnormal, add them to:\n")
cat(delete_list_csv, "\n")

cat("\nExpected files_to_delete.csv format:\n")
cat("filename_with_extension,reason\n")
cat("pnt1_0005.sed,bad NIR\n")

cat("\nReading files_to_delete.csv and excluding listed files from export.\n")

# Read the user-edited CSV that lists files to exclude from export.
# The path to this CSV is defined above as delete_list_csv.
files_to_delete <- read_csv(
  delete_list_csv,
  show_col_types = FALSE
)

# Confirm that files_to_delete.csv has the expected columns.
if (!all(c("filename_with_extension", "reason") %in% names(files_to_delete))) {
  stop(
    "files_to_delete.csv must contain columns: filename_with_extension, reason"
  )
}

# Print the filenames that will be excluded from export.
cat("\nFiles listed for exclusion:\n")

if (nrow(files_to_delete) == 0) {
  cat("None\n")
} else {
  files_to_delete %>%
    select(filename_with_extension, reason) %>%
    print(n = Inf)
}

# Check whether any filenames listed for exclusion do not match
# files that were parsed successfully from the raw data folder.
unknown_delete_files <- files_to_delete %>%
  anti_join(
    valid_files %>% select(filename_with_extension = file),
    by = "filename_with_extension"
  )

# Warn the user if files_to_delete.csv includes filenames that were not
# found among the valid parsed files. This usually means there is a typo,
# a missing extension, or the file did not pass filename parsing.
if (nrow(unknown_delete_files) > 0) {
  warning("Some files listed in files_to_delete.csv were not found among valid files.")
  print(unknown_delete_files, n = Inf)
}

# ============================================================
# Exclude user-listed bad files and export good files
# ============================================================
#
# This section does not delete original files.
# It copies good files into:
#
#   data-collection-QC/outputs/good_files_full_filenames/
#
# with full harmonization filenames:
#
#   PIdataHarmonization2026_<herbariumCode>_<targetID>_<kitNumber>_<IDX>.<ext>

good_files <- valid_files %>%
  anti_join(
    files_to_delete %>% select(filename_with_extension),
    by = c("file" = "filename_with_extension")
  ) %>%
  mutate(
    export_projectID = "PIdataHarmonization2026",
    export_herbariumCode = herbariumCode,
    export_kitNumber = kitNumber_parsed_or_entered,
    
    export_file = paste0(
      export_projectID, "_",
      export_herbariumCode, "_",
      targetID, "_",
      export_kitNumber, "_",
      IDX, ".",
      extension
    ),
    
    export_path = file.path(good_files_dir, export_file)
  )

duplicate_export_files <- good_files %>%
  count(export_file, name = "n") %>%
  filter(n > 1)

if (nrow(duplicate_export_files) > 0) {
  print(duplicate_export_files, n = Inf)
  stop(
    "Duplicate export filenames detected. Check targetID, kitNumber, and IDX values before exporting."
  )
}

copy_success <- file.copy(
  from = good_files$path,
  to = good_files$export_path,
  overwrite = TRUE
)

good_files_export_log <- good_files %>%
  mutate(copy_success = copy_success) %>%
  select(
    filename_with_extension = file,
    export_file,
    targetID,
    kitNumber_parsed_or_entered,
    IDX,
    extension,
    path,
    export_path,
    copy_success
  )

write_csv(
  good_files_export_log,
  file.path(out_dir, "good_files_export_log.csv")
)

cat("\nGood files exported to:\n")
cat(good_files_dir, "\n")

cat("\nNumber of files listed for exclusion:", nrow(files_to_delete), "\n")
cat("Number of good files exported:", sum(copy_success), "\n")

if (any(!copy_success)) {
  warning("Some files failed to copy. See good_files_export_log.csv")
}

# ============================================================
# Done
# ============================================================

cat("\nQC complete.\n")
cat("Outputs written to:\n")
cat(out_dir, "\n")

cat("\nKey output files:\n")
cat(file.path(out_dir, "parsed_file_inventory.csv"), "\n")
cat(file.path(out_dir, "file_counts_by_material.csv"), "\n")
cat(file.path(out_dir, "bad_filenames.csv"), "\n")
cat(file.path(out_dir, "files_to_delete.csv"), "\n")
cat(file.path(out_dir, "good_files_export_log.csv"), "\n")
cat(all_qc_pdf, "\n")
cat(good_files_dir, "\n")
