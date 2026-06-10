# spectra-harmonization

A collaborative repository for the IHerbSpec spectral data harmonization project. 

Currently, the repo hosts the scripts to run the shiny app for spectral data quality control (QC). See instructions and functionality below.

Additional subdirectories for data transformation and analysis can be added to the `analyses/` directory as the project develops.

## Repository structure

```text
spectra-harmonization/
├── README.md
├── spectra-harmonization.Rproj
├── spectra-files-qc-app.R   # Shiny app for interactive QC and metadata export
├── scripts/
│   ├── qc_spectra_files.R   # Batch QC script (plots and file counts only)
│   └── qc_functions.R       # Shared helper functions (sourced automatically)
└── analyses/                # Open folder for shared analysis workflows
```

## Running the 

Open the repository in RStudio by double-clicking `spectra-harmonization.Rproj`. This sets the working directory to the repository root, which all scripts assume.

Install required packages if needed:

```r
install.packages(c("shiny", "shinyFiles", "tidyverse", "DT", "plotly", "fs","spectrolab","ggrepel"))
```

## Filename format

The scripts will read spectral files named using either the simple filename convention or the full file convention below. The material identifier must be one of the 19 valid materials listed below.

The script accepts SVC `.sig`, SE `.sed`, ASD `.asd`, and `.txt` files, which are read with the spectrolab package.

### Simple local filename
```
[TC]<material>_<IDX>.<ext>
```
Examples: `fab2_0000.sig`, `TCfab2_0000.sed`

The `TC` prefix on the material is optional and will be stripped internally.

### Full harmonization filename
```
PIdataHarmonization2026_HC<herbariumCode>_TC<material>_kit<n>_<IDX>.<ext>
```
Example: `PIdataHarmonization2026_HCHUH_TCfab2_kit1_0000.sig`

### Valid materials

All files need to contain these exact materials IDs or they will be flagged as unrecognized.

| | | | |
|---|---|---|---|
| fab2 | fab5 | fel2 | fel3 |
| magmac-ab | magmac-ad | pap11 | pap6 |
| phymac-ab | phymac-ad | pnt1 | pnt2 |
| pnt3 | pnt4 | ravmad-ab | ravmad-ad |
| tcb | tcw | tvk | |


---

## Shiny QC and metadata app (`spectra-files-qc-app.R`)

The Shiny app supports interactive visual QC, file flagging, metadata capture, and standardized file export.

### Running the app

```r
shiny::runApp("spectra-files-qc-app.R")
```

### Workflow

1. **Select a directory** containing your raw spectral files.
2. **Review file checks** — the app reports unrecognised filenames, materials with fewer than five measurements, and any of the 19 expected materials not present in the file set.
3. **Review QC plots** — spectra are plotted per material. Click a spectrum line or use the checkboxes to flag individual files for removal.
4. **Fill in the metadata form** — instrument, optical setup, operator, and measurement settings.
5. **Export** — unflagged files are copied to `shiny_qc_outputs/good_files_full_filenames/` using the full harmonization filename format. A metadata CSV is written alongside.

### Outputs

All outputs are written to `shiny_qc_outputs/` inside the selected spectral files directory (not the repository root directory):

| File | Description |
|------|-------------|
| `good_files_full_filenames/` | Renamed good files in full harmonization format |
| `metadata-<projectId>_<herbariumCode>.csv` | Per-file metadata with a blank `comment` column |
| `good_files_export_log.csv` | Mapping from original to exported filename |
| `file_counts_by_material.csv` | File count per material |
| `bad_filenames.csv` | Files that failed filename validation |
| `missing_materials_warning.csv` | Materials from the expected 19 not found (if any) |
| `files_flagged_for_removal.csv` | Files flagged during visual QC |
| `qc_plots/` | Individual PNG QC plots per material |
| `ALL_materials_QC_labeled_spectra.pdf` | Multipage QC plot PDF |

---

## Batch QC script (`qc_spectra_files.R`)

Optionally, this is a non-interactive alternative for generating QC plots and file counts.

Create a `raw_data_files/` folder and an `outputs/` folder at the repository root. Place the spectral files within the `raw_data_files/` folder, then run the script.

The script will:
- Validate filenames and report unrecognised files to `bad_filenames.csv`
- Count files per material and flag any below the minimum threshold
- Report any of the 19 expected materials not found in the file set
- Generate per-material QC plots and a multipage PDF

---

## analyses/

An open folder for future work where participating researchers to share data transformation scripts, exploratory analyses, visualisations, and downstream harmonization workflows.

---

## Notes

Raw spectral files and outputs should generally not be committed to the repository. Add `raw_data_files/` and `outputs/` to `.gitignore` if working locally.
