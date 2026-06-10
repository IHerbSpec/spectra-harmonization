# ============================================================
# qc_functions.R
#
# Helper functions for qc_spectra_files.R.
# Sourced automatically at the top of that script.
# ============================================================

VALID_MATERIALS <- c(
  "fab2", "fel3", "fab5", "fel2",
  "magmac-ab", "magmac-ad",
  "pap11", "pap6",
  "phymac-ab", "phymac-ad",
  "pnt1", "pnt2", "pnt3", "pnt4",
  "ravmad-ab", "ravmad-ad",
  "tcb", "tcw", "tvk"
)


# ---- List spectral files -----------------------------------------------
#
# Returns a tibble with one row per file found in spectra_dir.
# Stops if no files with valid extensions are found.

list_spectral_files <- function(spectra_dir) {
  valid_extensions <- c("sed", "sig", "asd", "txt")

  spectral_files <- list.files(
    spectra_dir,
    pattern     = paste0("\\.(", paste(valid_extensions, collapse = "|"), ")$"),
    full.names  = TRUE,
    ignore.case = TRUE
  )

  if (length(spectral_files) == 0) {
    stop(
      "No spectral files found in: ", spectra_dir,
      "\nExpected extensions: .sed, .sig, .asd, .txt"
    )
  }

  tibble(
    path            = spectral_files,
    file            = basename(spectral_files),
    extension       = str_to_lower(str_extract(basename(spectral_files), "[^.]+$")),
    filename_no_ext = str_remove(basename(spectral_files), "\\.[^.]+$")
  )
}


# ---- Parse filenames ---------------------------------------------------
#
# Applies regex matching for two supported filename structures:
#
#   Full harmonization: PIdataHarmonization2026_HC<code>_TC<material>_kit<n>_IDX.ext
#   Simple local:       [TC]<material>_IDX.ext
#
# The optional TC prefix on the material is stripped before validation and
# analysis; the canonical targetID is always the bare material string (e.g.
# "fab2", not "TCfab2"). targetID must be one of the 19 VALID_MATERIALS.
#
# Returns a list:
#   $valid  — files that matched, parsed cleanly, and have a valid material
#   $bad    — files that did not match, are missing required fields, or have
#             an unrecognised material

parse_filenames <- function(file_inventory, kitNumber = NA_character_) {
  shared_regex <- "^(PIdataHarmonization2026)_([A-Za-z0-9]+)_(.+)_(\\d+)_(\\d+)$"
  local_regex  <- "^(.+)_(\\d+)$"

  parsed <- file_inventory |>
    mutate(
      is_shared_name = str_detect(filename_no_ext, shared_regex),
      is_local_name  = str_detect(filename_no_ext, local_regex) & !is_shared_name,
      filename_ok    = is_shared_name | is_local_name
    ) |>
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
      targetID_local  = if_else(
        is_local_name,
        str_match(filename_no_ext, local_regex)[, 2],
        NA_character_
      ),
      kitNumber_local = if_else(is_local_name, kitNumber, NA_character_),
      IDX_local       = if_else(
        is_local_name,
        str_match(filename_no_ext, local_regex)[, 3],
        NA_character_
      )
    ) |>
    mutate(
      targetID_raw                = coalesce(targetID_shared, targetID_local),
      targetID                    = str_remove(targetID_raw, "^TC"),
      kitNumber_parsed_or_entered = coalesce(kitNumber_shared, kitNumber_local),
      IDX                         = coalesce(IDX_shared, IDX_local),
      kitNumber_int               = as.integer(kitNumber_parsed_or_entered),
      IDX_int                     = as.integer(IDX),
      material_valid              = targetID %in% VALID_MATERIALS,
      material_id                 = paste(targetID, kitNumber_parsed_or_entered, sep = "_"),
      filename_type               = case_when(
        is_shared_name ~ "full_harmonization_filename",
        is_local_name  ~ "simple_local_filename",
        TRUE           ~ "unrecognized_filename"
      )
    )

  list(
    valid = parsed |> filter(filename_ok, material_valid, !is.na(targetID), !is.na(kitNumber_int), !is.na(IDX_int)),
    bad   = parsed |> filter(!filename_ok | !material_valid | is.na(targetID) | is.na(kitNumber_int) | is.na(IDX_int))
  )
}


# ---- Count files per material ------------------------------------------
#
# Returns a tibble with one row per targetID + kitNumber combination,
# showing file counts and a QC flag for materials below min_expected_files.

check_file_counts <- function(valid_files, min_expected_files) {
  valid_files |>
    count(targetID, name = "n_files") |>
    mutate(
      qc_flag = if_else(
        n_files < min_expected_files,
        paste0("FLAG: fewer than ", min_expected_files, " files"),
        "OK"
      )
    ) |>
    arrange(qc_flag, targetID)
}

check_missing_materials <- function(valid_files) {
  found <- unique(valid_files$targetID)
  tibble(missing_material = VALID_MATERIALS[!VALID_MATERIALS %in% found])
}


# ---- Read spectral files -----------------------------------------------
#
# Reads all valid spectral files by extension and combines them into a
# single spectrolab spectra object.

read_spectral_files <- function(valid_files) {
  valid_extensions <- c("sed", "sig", "asd", "txt")

  read_one_extension <- function(ext) {
    these_files <- valid_files |> filter(extension == ext)
    if (nrow(these_files) == 0) return(NULL)

    cat("Reading", nrow(these_files), paste0(".", ext), "files\n")
    read_spectra(
      path             = these_files$path,
      format           = ext,
      type             = "target_reflectance",
      extract_metadata = FALSE
    )
  }

  spec_list        <- map(valid_extensions, read_one_extension)
  names(spec_list) <- valid_extensions
  spec_list        <- compact(spec_list)

  if (length(spec_list) == 0) stop("No readable spectral files found.")

  spec <- if (length(spec_list) == 1) {
    spec_list[[1]]
  } else {
    reduce(spec_list, spectrolab::combine)
  }

  if (is.list(spec) && !inherits(spec, "spectra")) {
    stop(
      "read_spectra() returned a list, likely because files have incompatible wavelength grids. ",
      "Inspect or resample files before continuing."
    )
  }

  spec
}


# ---- Convert spectra to long format ------------------------------------
#
# Converts a spectrolab spectra object to a long tibble joined to
# parsed filename metadata from valid_files.

spectra_to_long <- function(spec, valid_files) {
  spec_df           <- as.data.frame(spec, fix_names = "none", metadata = FALSE)
  names(spec_df)[1] <- "sample_name"

  spec_df |>
    pivot_longer(
      cols      = -sample_name,
      names_to  = "wavelength",
      values_to = "reflectance"
    ) |>
    mutate(
      wavelength = as.numeric(wavelength),
      file       = basename(sample_name)
    ) |>
    left_join(
      valid_files |>
        select(
          file, extension, filename_type,
          projectID_parsed, herbariumCode_parsed,
          targetID, kitNumber_parsed_or_entered,
          kitNumber_int, IDX, IDX_int, material_id
        ),
      by = "file"
    ) |>
    filter(!is.na(material_id))
}


# ---- Plot one material --------------------------------------------------
#
# Returns a ggplot QC plot for one material_id.
# Each line is one file; labels show the filename without extension.

plot_material_qc <- function(spectra_long, material_name) {
  dat <- spectra_long |>
    filter(material_id == material_name) |>
    mutate(file_label = str_remove(file, paste0("\\.", extension, "$")))

  if (nrow(dat) == 0) return(NULL)

  label_df <- dat |>
    group_by(file, file_label) |>
    filter(wavelength == max(wavelength, na.rm = TRUE)) |>
    ungroup()

  ggplot(dat, aes(x = wavelength, y = reflectance, group = file, color = file)) +
    geom_line(linewidth = 0.55, alpha = 0.85, show.legend = FALSE) +
    ggrepel::geom_text_repel(
      data          = label_df,
      aes(label     = file_label),
      hjust         = 0,
      nudge_x       = 30,
      direction     = "y",
      size          = 3.2,
      segment.size  = 0.25,
      segment.alpha = 0.5,
      show.legend   = FALSE
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.25))) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title    = paste("QC plot:", material_name),
      subtitle = "Each line is one spectral file; labels identify filenames for visual QC",
      x        = "Wavelength (nm)",
      y        = "Reflectance"
    ) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
}


# ---- Save QC plots (PNGs + multipage PDF) ------------------------------
#
# Saves one PNG per material to qc_plot_dir, and a multipage PDF of all
# materials to out_dir. Returns the path to the PDF invisibly.

save_qc_plots <- function(spectra_long, qc_plot_dir, out_dir) {
  materials <- spectra_long |>
    distinct(material_id) |>
    arrange(material_id) |>
    pull(material_id)

  for (mat in materials) {
    p        <- plot_material_qc(spectra_long, mat)
    safe_mat <- str_replace_all(mat, "[^A-Za-z0-9]+", "_")
    ggsave(
      filename = file.path(qc_plot_dir, paste0("QC_", safe_mat, ".png")),
      plot     = p,
      width    = 10,
      height   = 6,
      dpi      = 300
    )
  }

  pdf_path <- file.path(out_dir, "ALL_materials_QC_labeled_spectra.pdf")
  pdf(pdf_path, width = 10, height = 6)
  for (mat in materials) print(plot_material_qc(spectra_long, mat))
  dev.off()

  invisible(pdf_path)
}


# ---- Export good files -------------------------------------------------
#
# Copies files that are not listed in files_to_delete into good_files_dir
# using the full harmonization filename convention:
#
#   PIdataHarmonization2026_HC<herbariumCode>_TC<targetID>_kit<kitNumber>_<IDX>.<ext>
#
# Writes good_files_export_log.csv to out_dir.
# Returns the export log tibble invisibly.

export_good_files <- function(valid_files, files_to_delete, herbariumCode,
                              good_files_dir, out_dir) {
  good_files <- valid_files |>
    anti_join(
      files_to_delete |> select(filename_with_extension),
      by = c("file" = "filename_with_extension")
    ) |>
    mutate(
      export_file = paste0(
        "PIdataHarmonization2026", "_",
        "HC", herbariumCode, "_",
        "TC", targetID, "_",
        "kit", kitNumber_parsed_or_entered, "_",
        IDX, ".",
        extension
      ),
      export_path = file.path(good_files_dir, export_file)
    )

  duplicate_export_files <- good_files |>
    count(export_file, name = "n") |>
    filter(n > 1)

  if (nrow(duplicate_export_files) > 0) {
    print(duplicate_export_files, n = Inf)
    stop(
      "Duplicate export filenames detected. ",
      "Check targetID, kitNumber, and IDX values before exporting."
    )
  }

  copy_success <- file.copy(
    from      = good_files$path,
    to        = good_files$export_path,
    overwrite = TRUE
  )

  export_log <- good_files |>
    mutate(copy_success = copy_success) |>
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

  write_csv(export_log, file.path(out_dir, "good_files_export_log.csv"))

  if (any(!copy_success)) warning("Some files failed to copy. See good_files_export_log.csv")

  invisible(export_log)
}
