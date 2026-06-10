# ============================================================
# IHerbSpec data harmonization QC + metadata Shiny app
#
# Place this file at the root of the spectra-harmonization repo
# and run with:
#   shiny::runApp("app.R")
#
# Required packages:
#   install.packages(c("shiny", "shinyFiles", "tidyverse", "DT", "plotly", "fs"))
#   remotes::install_github("plant-data/spectrolab") # if spectrolab is not installed
# ============================================================

library(shiny)
library(shinyFiles)
library(tidyverse)
library(DT)
library(plotly)
library(fs)
library(spectrolab)

# ============================================================
# Constants
# ============================================================

VALID_EXTENSIONS <- c("sed", "sig", "asd", "txt")
MIN_EXPECTED_FILES <- 5

DEFAULT_PROJECT_ID <- "dataHarmonization2026"
DEFAULT_BACKGROUND_CLASS <- "BGB"
DEFAULT_HAS_LOW_REFLECTANCE_BACKGROUND <- TRUE
DEFAULT_HAS_BACKGROUND_IN_MEASUREMENT <- FALSE
DEFAULT_PERCENT_BACKGROUND_IN_MEASUREMENT <- 0
DEFAULT_BACKGROUND_DESCRIPTION <- "Musou IR Flock"

VALID_MATERIALS <- c(
  "fab2", "fel3", "fab5", "fel2",
  "magmac-ab", "magmac-ad",
  "pap11", "pap6",
  "phymac-ab", "phymac-ad",
  "pnt1", "pnt2", "pnt3", "pnt4",
  "ravmad-ab", "ravmad-ad",
  "tcb", "tcw", "tvk"
)

# ============================================================
# Helper functions adapted from data-collection-QC/scripts/qc_functions.R
# ============================================================

list_spectral_files <- function(spectra_dir) {
  spectral_files <- list.files(
    spectra_dir,
    pattern = paste0("\\.(", paste(VALID_EXTENSIONS, collapse = "|"), ")$"),
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(spectral_files) == 0) {
    stop(
      "No spectral files found in: ", spectra_dir,
      "\nExpected extensions: .", paste(VALID_EXTENSIONS, collapse = ", .")
    )
  }

  tibble(
    path = spectral_files,
    short_filename = basename(spectral_files),
    extension = str_to_lower(str_extract(basename(spectral_files), "[^.]+$")),
    filename_no_ext = str_remove(basename(spectral_files), "\\.[^.]+$")
  )
}

parse_filenames <- function(file_inventory, kitNumber_default = NA_character_) {
  # New full harmonization filename:
  #   projectId_herbariumCode_targetId_kitNumber_measurementIndex.ext
  # Simple local filename:
  #   targetId_measurementIndex.ext
  #
  # This parser allows underscores inside targetId by anchoring the final two
  # fields for the full format and the final field for the local format.

  full_regex <- "^([^_]+)_([A-Za-z0-9]+)_(.+)_([0-9]+)_([0-9]+)$"
  local_regex <- "^(.+)_([0-9]+)$"

  parsed <- file_inventory |>
    mutate(
      is_full_name = str_detect(filename_no_ext, full_regex),
      is_local_name = str_detect(filename_no_ext, local_regex) & !is_full_name,
      filename_ok = is_full_name | is_local_name,
      projectId_parsed = if_else(
        is_full_name,
        str_match(filename_no_ext, full_regex)[, 2],
        NA_character_
      ),
      herbariumCode_parsed = if_else(
        is_full_name,
        str_match(filename_no_ext, full_regex)[, 3],
        NA_character_
      ),
      targetId_full = if_else(
        is_full_name,
        str_match(filename_no_ext, full_regex)[, 4],
        NA_character_
      ),
      kitNumber_full = if_else(
        is_full_name,
        str_match(filename_no_ext, full_regex)[, 5],
        NA_character_
      ),
      measurementIndex_full = if_else(
        is_full_name,
        str_match(filename_no_ext, full_regex)[, 6],
        NA_character_
      ),
      targetId_local = if_else(
        is_local_name,
        str_match(filename_no_ext, local_regex)[, 2],
        NA_character_
      ),
      kitNumber_local = if_else(is_local_name, as.character(kitNumber_default), NA_character_),
      measurementIndex_local = if_else(
        is_local_name,
        str_match(filename_no_ext, local_regex)[, 3],
        NA_character_
      ),
      targetId_raw = coalesce(targetId_full, targetId_local),
      targetId = str_remove(targetId_raw, "^TC"),
      kitNumber = coalesce(kitNumber_full, kitNumber_local),
      measurementIndex = coalesce(measurementIndex_full, measurementIndex_local),
      kitNumber_int = suppressWarnings(as.integer(kitNumber)),
      measurementIndex_int = suppressWarnings(as.integer(measurementIndex)),
      material_valid = targetId %in% VALID_MATERIALS,
      material_id = paste(targetId, kitNumber, sep = "_"),
      filename_type = case_when(
        is_full_name ~ "full_harmonization_filename",
        is_local_name ~ "simple_local_filename",
        TRUE ~ "unrecognized_filename"
      )
    )

  list(
    valid = parsed |>
      filter(filename_ok, material_valid, !is.na(targetId), !is.na(kitNumber_int), !is.na(measurementIndex_int)),
    bad = parsed |>
      filter(!filename_ok | !material_valid | is.na(targetId) | is.na(kitNumber_int) | is.na(measurementIndex_int))
  )
}

check_file_counts <- function(valid_files, min_expected_files = MIN_EXPECTED_FILES) {
  valid_files |>
    count(targetId, kitNumber, material_id, name = "n_files") |>
    mutate(
      qc_flag = if_else(
        n_files < min_expected_files,
        paste0("FLAG: fewer than ", min_expected_files, " files"),
        "OK"
      )
    ) |>
    arrange(desc(qc_flag != "OK"), targetId, kitNumber)
}

check_missing_materials <- function(valid_files) {
  found <- unique(valid_files$targetId)
  tibble(missing_material = VALID_MATERIALS[!VALID_MATERIALS %in% found])
}

read_spectral_files <- function(valid_files) {
  read_one_extension <- function(ext) {
    these_files <- valid_files |> filter(extension == ext)
    if (nrow(these_files) == 0) return(NULL)

    read_spectra(
      path = these_files$path,
      format = ext,
      type = "target_reflectance",
      extract_metadata = FALSE
    )
  }

  spec_list <- map(VALID_EXTENSIONS, read_one_extension)
  spec_list <- purrr::compact(spec_list)

  if (length(spec_list) == 0) stop("No readable spectral files found.")

  spec <- if (length(spec_list) == 1) {
    spec_list[[1]]
  } else {
    purrr::reduce(spec_list, spectrolab::combine)
  }

  if (is.list(spec) && !inherits(spec, "spectra")) {
    stop(
      "read_spectra() returned a list, likely because files have incompatible wavelength grids. ",
      "Inspect or resample files before continuing."
    )
  }

  spec
}

spectra_to_long <- function(spec, valid_files) {
  spec_df <- as.data.frame(spec, fix_names = "none", metadata = FALSE)
  names(spec_df)[1] <- "sample_name"

  spec_df |>
    pivot_longer(
      cols = -sample_name,
      names_to = "wavelength",
      values_to = "reflectance"
    ) |>
    mutate(
      wavelength = as.numeric(wavelength),
      short_filename = basename(sample_name)
    ) |>
    left_join(
      valid_files |>
        select(
          short_filename, extension, filename_type,
          projectId_parsed, herbariumCode_parsed,
          targetId, kitNumber, kitNumber_int,
          measurementIndex, measurementIndex_int,
          material_id
        ),
      by = "short_filename"
    ) |>
    filter(!is.na(material_id))
}

plot_material_static <- function(spectra_long, material_name) {
  dat <- spectra_long |>
    filter(material_id == material_name) |>
    mutate(file_label = str_remove(short_filename, paste0("\\.", extension, "$")))

  ggplot(dat, aes(x = wavelength, y = reflectance, group = short_filename, color = short_filename)) +
    geom_line(linewidth = 0.55, alpha = 0.85, show.legend = FALSE) +
    labs(
      title = paste("QC plot:", material_name),
      subtitle = "Each line is one spectral file",
      x = "Wavelength (nm)",
      y = "Reflectance"
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
}

plot_material_interactive <- function(spectra_long, material_name, flagged_files = character()) {
  dat <- spectra_long |>
    filter(material_id == material_name) |>
    arrange(short_filename, wavelength)

  p <- plot_ly(source = "qcplot")

  for (fn in unique(dat$short_filename)) {
    one <- dat |> filter(short_filename == fn)
    line_width <- if (fn %in% flagged_files) 4 else 1.5
    p <- add_trace(
      p,
      data = one,
      x = ~wavelength,
      y = ~reflectance,
      type = "scatter",
      mode = "lines",
      name = fn,
      customdata = ~short_filename,
      hovertemplate = paste0(
        "File: %{customdata}<br>",
        "Wavelength: %{x}<br>",
        "Reflectance: %{y}<extra></extra>"
      ),
      line = list(width = line_width)
    )
  }

  p |>
    layout(
      title = paste("QC plot:", material_name),
      xaxis = list(title = "Wavelength (nm)"),
      yaxis = list(title = "Reflectance", range = c(0, 1)),
      legend = list(orientation = "v")
    ) |>
    event_register("plotly_click")
}

save_qc_plots <- function(spectra_long, qc_plot_dir, out_dir) {
  dir_create(qc_plot_dir)
  dir_create(out_dir)

  materials <- spectra_long |>
    distinct(material_id) |>
    arrange(material_id) |>
    pull(material_id)

  for (mat in materials) {
    p <- plot_material_static(spectra_long, mat)
    safe_mat <- str_replace_all(mat, "[^A-Za-z0-9]+", "_")
    ggsave(
      filename = file.path(qc_plot_dir, paste0("QC_", safe_mat, ".png")),
      plot = p,
      width = 10,
      height = 6,
      dpi = 300
    )
  }

  pdf_path <- file.path(out_dir, "ALL_materials_QC_labeled_spectra.pdf")
  pdf(pdf_path, width = 10, height = 6)
  on.exit(dev.off(), add = TRUE)
  for (mat in materials) print(plot_material_static(spectra_long, mat))

  pdf_path
}

infer_measurement_settings <- function(instrumentModel) {
  model <- str_to_lower(instrumentModel %||% "")
  case_when(
    str_detect(model, "hr-?1024i|svc") ~ "Low light, 3 second integration time",
    str_detect(model, "naturaspec|natura spec") ~ "Automatic integration, 40 averaging time.",
    TRUE ~ NA_character_
  )
}

build_full_filename_no_ext <- function(projectId, herbariumCode, targetId, kitNumber, measurementIndex) {
  paste(paste0("PI", projectId), paste0("HC", herbariumCode), paste0("TC", targetId), paste0("kit", kitNumber), measurementIndex, sep = "_")
}

export_good_files_and_metadata <- function(valid_files, flagged_files, form_values, out_dir, good_files_dir) {
  dir_create(out_dir)
  dir_create(good_files_dir)

  good_files <- valid_files |>
    filter(!short_filename %in% flagged_files) |>
    mutate(
      projectId = form_values$projectId,
      herbariumCode = form_values$herbariumCode,
      kitNumber = as.character(kitNumber),
      measurementIndex = as.character(measurementIndex),
      full_filename = build_full_filename_no_ext(
        projectId = projectId,
        herbariumCode = herbariumCode,
        targetId = targetId,
        kitNumber = kitNumber,
        measurementIndex = measurementIndex
      ),
      export_file = paste0(full_filename, ".", extension),
      export_path = file.path(good_files_dir, export_file)
    )

  duplicate_export_files <- good_files |>
    count(export_file, name = "n") |>
    filter(n > 1)

  if (nrow(duplicate_export_files) > 0) {
    stop("Duplicate export filenames detected. Check targetId, kitNumber, and measurementIndex values.")
  }

  copy_success <- file.copy(
    from = good_files$path,
    to = good_files$export_path,
    overwrite = TRUE
  )

  export_log <- good_files |>
    mutate(copy_success = copy_success) |>
    select(short_filename, export_file, targetId, kitNumber, measurementIndex, path, export_path, copy_success)

  write_csv(export_log, file.path(out_dir, "good_files_export_log.csv"))

  metadata <- good_files |>
    transmute(
      full_filename,
      short_filename,
      targetId,
      measurementIndex,
      projectId = form_values$projectId,
      herbariumCode = form_values$herbariumCode,
      kitNumber = as.numeric(kitNumber),
      backgroundClass = form_values$backgroundClass,
      hasLowReflectanceBackground = form_values$hasLowReflectanceBackground,
      hasBackgroundInMeasurement = form_values$hasBackgroundInMeasurement,
      percentBackgroundInMeasurement = form_values$percentBackgroundInMeasurement,
      backgroundDescription = form_values$backgroundDescription,
      whiteReferenceDescription = form_values$whiteReferenceDescription,
      instrumentModel = form_values$instrumentModel,
      opticalSetupDescription = form_values$opticalSetupDescription,
      measurementSettings = form_values$measurementSettings,
      operator = form_values$operator,
      instrumentCalibrationDate = form_values$instrumentCalibrationDate,
      lightSourceType = form_values$lightSourceType,
      distanceTargetToSensor = form_values$distanceTargetToSensor,
      lensFieldOfView = form_values$lensFieldOfView,
      angleLightToSensor = form_values$angleLightToSensor,
      measurementAreaDiameter = form_values$measurementAreaDiameter,
      comment = ""
    ) |>
    arrange(targetId, as.integer(measurementIndex))

  metadata_path <- file.path(out_dir, paste0("metadata-", form_values$projectId, "_", form_values$herbariumCode, ".csv"))
  write_csv(metadata, metadata_path)

  files_to_delete <- tibble(
    filename_with_extension = flagged_files,
    reason = "flagged in Shiny app"
  )
  write_csv(files_to_delete, file.path(out_dir, "files_flagged_for_removal.csv"))

  list(
    good_files = good_files,
    export_log = export_log,
    metadata = metadata,
    metadata_path = metadata_path,
    flagged_path = file.path(out_dir, "files_flagged_for_removal.csv"),
    good_files_dir = good_files_dir
  )
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

# ============================================================
# Shiny UI
# ============================================================

ui <- fluidPage(
  titlePanel("IHerbSpec spectra QC and metadata builder"),

  sidebarLayout(
    sidebarPanel(
      h4("1. Select spectral-file directory"),
      shinyDirButton("choose_dir", "Choose directory", "Select the folder containing .sed, .sig, .asd, or .txt files"),
      verbatimTextOutput("selected_dir"),
      numericInput("kitNumber", "Kit number used for simple local filenames", value = 1, min = 1, step = 1),
      actionButton("read_files", "Read files and generate QC plots", class = "btn-primary"),
      hr(),
      h4("2. Required metadata"),
      textInput("herbariumCode", "herbariumCode", value = "HUH"),
      textInput("whiteReferenceDescription", "whiteReferenceDescription", value = "Spectralon"),
      textInput("instrumentModel", "instrumentModel", value = ""),
      textInput("opticalSetupDescription", "opticalSetupDescription", value = ""),
      uiOutput("measurement_settings_ui"),
      hr(),
      h4("3. Confirm auto-populated fields"),
      textInput("projectId", "projectId", value = DEFAULT_PROJECT_ID),
      textInput("backgroundClass", "backgroundClass", value = DEFAULT_BACKGROUND_CLASS),
      checkboxInput("hasLowReflectanceBackground", "hasLowReflectanceBackground", value = DEFAULT_HAS_LOW_REFLECTANCE_BACKGROUND),
      checkboxInput("hasBackgroundInMeasurement", "hasBackgroundInMeasurement", value = DEFAULT_HAS_BACKGROUND_IN_MEASUREMENT),
      numericInput("percentBackgroundInMeasurement", "percentBackgroundInMeasurement", value = DEFAULT_PERCENT_BACKGROUND_IN_MEASUREMENT, min = 0, max = 100),
      textInput("backgroundDescription", "backgroundDescription", value = DEFAULT_BACKGROUND_DESCRIPTION),
      hr(),
      h4("4. Optional metadata"),
      textInput("operator", "operator", value = ""),
      textInput("instrumentCalibrationDate", "instrumentCalibrationDate", value = ""),
      textInput("lightSourceType", "lightSourceType", value = ""),
      textInput("distanceTargetToSensor", "distanceTargetToSensor", value = ""),
      textInput("lensFieldOfView", "lensFieldOfView", value = ""),
      textInput("angleLightToSensor", "angleLightToSensor", value = ""),
      textInput("measurementAreaDiameter", "measurementAreaDiameter", value = ""),
      textAreaInput("globalComment", "Optional note for this session; not copied to row-level comment", value = "", rows = 2),
      hr(),
      actionButton("export", "Export good files and metadata CSV", class = "btn-success")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel(
          "File checks",
          h4("Missing materials (not found in any filename)"),
          DTOutput("missing_materials"),
          h4("Materials with fewer than five samples"),
          DTOutput("flagged_materials"),
          h4("All material counts"),
          DTOutput("file_counts"),
          h4("Bad or unrecognized filenames"),
          DTOutput("bad_filenames"),
          h4("QC output status"),
          verbatimTextOutput("qc_status")
        ),
        tabPanel(
          "Visual QC",
          fluidRow(
            column(5, selectInput("material", "Material", choices = character())),
            column(7, checkboxGroupInput("flagged_for_material", "Flag files for removal", choices = character()))
          ),
          plotlyOutput("qc_plot", height = "650px"),
          helpText("Click a spectrum line to toggle that file as flagged for removal. You can also use the checkboxes above."),
          h4("All files flagged for removal"),
          DTOutput("flagged_files")
        ),
        tabPanel(
          "Metadata preview",
          p("The exported metadata CSV includes an additional blank comment field. Users can fill this out manually after export."),
          DTOutput("metadata_preview"),
          h4("Export status"),
          verbatimTextOutput("export_status")
        )
      )
    )
  )
)

# ============================================================
# Shiny server
# ============================================================

server <- function(input, output, session) {
  volumes <- c(Home = fs::path_home(), WorkingDirectory = getwd())
  shinyDirChoose(input, "choose_dir", roots = volumes, session = session, restrictions = system.file(package = "base"))

  spectra_dir <- reactiveVal(NULL)
  parsed_data <- reactiveVal(NULL)
  spectra_long_data <- reactiveVal(NULL)
  pdf_path_val <- reactiveVal(NULL)
  flagged_files <- reactiveVal(character())
  export_result <- reactiveVal(NULL)

  observeEvent(input$choose_dir, {
    selected <- parseDirPath(volumes, input$choose_dir)
    spectra_dir(selected)
  })

  output$selected_dir <- renderText({
    req(spectra_dir())
    spectra_dir()
  })

  output$measurement_settings_ui <- renderUI({
    inferred <- infer_measurement_settings(input$instrumentModel)
    if (is.na(inferred)) {
      textInput("measurementSettings", "measurementSettings", value = "")
    } else {
      textInput("measurementSettings", "measurementSettings", value = inferred)
    }
  })

  observeEvent(input$read_files, {
    req(spectra_dir())

    withProgress(message = "Reading spectral files and generating QC plots", value = 0, {
      incProgress(0.10, detail = "Listing and parsing filenames")
      file_inventory <- list_spectral_files(spectra_dir())
      parsed <- parse_filenames(file_inventory, kitNumber_default = as.character(input$kitNumber))

      out_dir <- file.path(spectra_dir(), "shiny_qc_outputs")
      qc_plot_dir <- file.path(out_dir, "qc_plots")
      good_files_dir <- file.path(out_dir, "good_files_full_filenames")
      dir_create(out_dir)
      dir_create(qc_plot_dir)
      dir_create(good_files_dir)

      write_csv(parsed$bad, file.path(out_dir, "bad_filenames.csv"))
      write_csv(parsed$valid, file.path(out_dir, "parsed_file_inventory.csv"))

      if (nrow(parsed$valid) == 0) stop("No files matched the expected filename conventions.")

      incProgress(0.35, detail = "Reading spectra")
      spec <- read_spectral_files(parsed$valid)
      spectra_long <- spectra_to_long(spec, parsed$valid)

      incProgress(0.35, detail = "Saving multipage PDF and individual plots")
      pdf_path <- save_qc_plots(spectra_long, qc_plot_dir, out_dir)

      parsed_data(list(
        valid = parsed$valid,
        bad = parsed$bad,
        file_counts = check_file_counts(parsed$valid),
        missing_materials = check_missing_materials(parsed$valid),
        out_dir = out_dir,
        qc_plot_dir = qc_plot_dir,
        good_files_dir = good_files_dir
      ))
      spectra_long_data(spectra_long)
      pdf_path_val(pdf_path)
      flagged_files(character())

      materials <- spectra_long |>
        distinct(material_id) |>
        arrange(material_id) |>
        pull(material_id)
      updateSelectInput(session, "material", choices = materials, selected = materials[1])

      incProgress(0.20, detail = "Done")
    })
  })

  output$file_counts <- renderDT({
    req(parsed_data())
    datatable(parsed_data()$file_counts, options = list(pageLength = 10))
  })

  output$missing_materials <- renderDT({
    req(parsed_data())
    datatable(parsed_data()$missing_materials, options = list(pageLength = 25))
  })

  output$flagged_materials <- renderDT({
    req(parsed_data())
    flagged <- parsed_data()$file_counts |> filter(n_files < MIN_EXPECTED_FILES)
    datatable(flagged, options = list(pageLength = 10))
  })

  output$bad_filenames <- renderDT({
    req(parsed_data())
    datatable(parsed_data()$bad, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$qc_status <- renderText({
    req(parsed_data(), pdf_path_val())
    paste(
      "Multipage PDF:", pdf_path_val(),
      "\nIndividual PNGs:", parsed_data()$qc_plot_dir,
      "\nParsed inventory:", file.path(parsed_data()$out_dir, "parsed_file_inventory.csv"),
      "\nFile counts:", file.path(parsed_data()$out_dir, "file_counts_by_material.csv"),
      sep = ""
    )
  })

  observeEvent(parsed_data(), {
    write_csv(parsed_data()$file_counts, file.path(parsed_data()$out_dir, "file_counts_by_material.csv"))
    if (nrow(parsed_data()$missing_materials) > 0) {
      write_csv(parsed_data()$missing_materials, file.path(parsed_data()$out_dir, "missing_materials_warning.csv"))
    }
  })

  observeEvent(input$material, {
    req(parsed_data(), input$material)
    choices <- parsed_data()$valid |>
      filter(material_id == input$material) |>
      arrange(short_filename) |>
      pull(short_filename)
    updateCheckboxGroupInput(
      session,
      "flagged_for_material",
      choices = choices,
      selected = intersect(choices, flagged_files())
    )
  }, ignoreInit = TRUE)

  observeEvent(input$flagged_for_material, {
    req(parsed_data(), input$material)
    choices <- parsed_data()$valid |>
      filter(material_id == input$material) |>
      pull(short_filename)

    old <- flagged_files()
    new <- union(setdiff(old, choices), input$flagged_for_material %||% character())
    flagged_files(sort(unique(new)))
  }, ignoreInit = TRUE)

  observeEvent(event_data("plotly_click", source = "qcplot"), {
    click <- event_data("plotly_click", source = "qcplot")
    req(click$customdata)
    fn <- click$customdata[[1]]
    old <- flagged_files()
    new <- if (fn %in% old) setdiff(old, fn) else union(old, fn)
    flagged_files(sort(unique(new)))

    if (!is.null(input$material) && !is.null(parsed_data())) {
      choices <- parsed_data()$valid |>
        filter(material_id == input$material) |>
        arrange(short_filename) |>
        pull(short_filename)
      updateCheckboxGroupInput(
        session,
        "flagged_for_material",
        choices = choices,
        selected = intersect(choices, flagged_files())
      )
    }
  }, ignoreInit = TRUE)

  output$qc_plot <- renderPlotly({
    req(spectra_long_data(), input$material)
    plot_material_interactive(spectra_long_data(), input$material, flagged_files())
  })

  output$flagged_files <- renderDT({
    tibble(filename_with_extension = flagged_files(), reason = "flagged in Shiny app") |>
      datatable(options = list(pageLength = 10))
  })

  current_form_values <- reactive({
    list(
      projectId = input$projectId,
      herbariumCode = input$herbariumCode,
      backgroundClass = input$backgroundClass,
      hasLowReflectanceBackground = input$hasLowReflectanceBackground,
      hasBackgroundInMeasurement = input$hasBackgroundInMeasurement,
      percentBackgroundInMeasurement = input$percentBackgroundInMeasurement,
      backgroundDescription = input$backgroundDescription,
      whiteReferenceDescription = input$whiteReferenceDescription,
      instrumentModel = input$instrumentModel,
      opticalSetupDescription = input$opticalSetupDescription,
      measurementSettings = input$measurementSettings,
      operator = input$operator,
      instrumentCalibrationDate = input$instrumentCalibrationDate,
      lightSourceType = input$lightSourceType,
      distanceTargetToSensor = input$distanceTargetToSensor,
      lensFieldOfView = input$lensFieldOfView,
      angleLightToSensor = input$angleLightToSensor,
      measurementAreaDiameter = input$measurementAreaDiameter,
      globalComment = input$globalComment
    )
  })

  metadata_preview_data <- reactive({
    req(parsed_data())
    form <- current_form_values()
    parsed_data()$valid |>
      filter(!short_filename %in% flagged_files()) |>
      transmute(
        full_filename = build_full_filename_no_ext(
          form$projectId,
          form$herbariumCode,
          targetId,
          kitNumber,
          measurementIndex
        ),
        short_filename,
        targetId,
        measurementIndex,
        projectId = form$projectId,
        herbariumCode = form$herbariumCode,
        kitNumber = as.numeric(kitNumber),
        backgroundClass = form$backgroundClass,
        hasLowReflectanceBackground = form$hasLowReflectanceBackground,
        hasBackgroundInMeasurement = form$hasBackgroundInMeasurement,
        percentBackgroundInMeasurement = form$percentBackgroundInMeasurement,
        backgroundDescription = form$backgroundDescription,
        whiteReferenceDescription = form$whiteReferenceDescription,
        instrumentModel = form$instrumentModel,
        opticalSetupDescription = form$opticalSetupDescription,
        measurementSettings = form$measurementSettings,
        operator = form$operator,
        instrumentCalibrationDate = form$instrumentCalibrationDate,
        lightSourceType = form$lightSourceType,
        distanceTargetToSensor = form$distanceTargetToSensor,
        lensFieldOfView = form$lensFieldOfView,
        angleLightToSensor = form$angleLightToSensor,
        measurementAreaDiameter = form$measurementAreaDiameter,
        comment = ""
      ) |>
      arrange(targetId, as.integer(measurementIndex))
  })

  output$metadata_preview <- renderDT({
    req(metadata_preview_data())
    datatable(metadata_preview_data(), options = list(pageLength = 10, scrollX = TRUE))
  })

  observeEvent(input$export, {
    req(parsed_data())

    required_missing <- c(
      herbariumCode = input$herbariumCode,
      kitNumber = as.character(input$kitNumber),
      whiteReferenceDescription = input$whiteReferenceDescription,
      instrumentModel = input$instrumentModel,
      opticalSetupDescription = input$opticalSetupDescription,
      measurementSettings = input$measurementSettings
    )
    required_missing <- names(required_missing)[is.na(required_missing) | required_missing == ""]

    if (length(required_missing) > 0) {
      showModal(modalDialog(
        title = "Required fields are missing",
        paste("Please fill:", paste(required_missing, collapse = ", ")),
        easyClose = TRUE
      ))
      return(NULL)
    }

    result <- export_good_files_and_metadata(
      valid_files = parsed_data()$valid,
      flagged_files = flagged_files(),
      form_values = current_form_values(),
      out_dir = parsed_data()$out_dir,
      good_files_dir = parsed_data()$good_files_dir
    )
    export_result(result)

    showModal(modalDialog(
      title = "Export complete",
      paste("Metadata CSV written to:", result$metadata_path),
      easyClose = TRUE
    ))
  })

  output$export_status <- renderText({
    req(export_result())
    result <- export_result()
    paste(
      "Metadata CSV: ", result$metadata_path,
      "\nGood files folder: ", result$good_files_dir,
      "\nFlagged-file CSV: ", result$flagged_path,
      "\nFiles exported: ", nrow(result$good_files),
      "\nFiles flagged for removal: ", length(flagged_files()),
      "\n\nThe metadata CSV includes a blank 'comment' field. Fill this out manually after export if needed.",
      sep = ""
    )
  })
}

shinyApp(ui, server)
