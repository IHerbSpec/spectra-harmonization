# spectra-harmonization

A collaborative repository for the IHerbSpec spectral data harmonization project. 

Currently, the repo hosts the scripts to run the shiny app for spectral data quality control (QC). See instructions and functionality below.

Additional subdirectories for data transformation and analysis can be added to the `analyses/` directory in the future as the project develops.

## Repository structure

```text
spectra-harmonization/
├── README.md
├── spectra-harmonization.Rproj
├── spectra-files-qc-app.R   # Shiny app for interactive QC and metadata export
├── test_data/               # .sed files for testing the Shiny app.
├── scripts/
│   ├── qc_spectra_files.R   # Batch QC script (plots and file counts only)
│   └── qc_functions.R       # Shared helper functions (sourced automatically)
└── analyses/                # Open folder for shared analysis workflows
```

## Running the Shiny spectral data QC and metadata builder app

This app was created to assist with project QC. It will: 
- read spectral files with simple filename convention
- create plots for visual inspection
- users can flag files for removal
- check the number of files for each kit material
- prompt input for required and optional metadata fields
- copy passing files to a new folder with converted full filename conventions
- generate an IHerbSpec-compatible metadata spreadsheet

### File and Filename format

The script accepts SVC `.sig`, SE `.sed`, ASD `.asd`, and `.txt` file extensions (`.<ext>`), which are read with the spectrolab package.

It is critical that each file contain the material identifier, which we treat in filenames as `targetClass`.

#### Valid material identifiers

All files need to contain these exact materials IDs bounded by underscores or the beginning of the filename, or they will be flagged as unrecognized.

| TC | TC | TC | TC | description|
|---|---|---|---|---|
| tcb | tcw | - | - | black background (provided) and white reference (not provided) |
| pnt1 | pnt2 | pnt3 | pnt4 | painted panels |
| fab2 | fab5 | - | - | fabrics |
| fel2 | fel3 | - | - | felts |
| tvk | - | - | - | tyvek |
| magmac-ab | magmac-ad | - | - | *Magnolia* leaf, abaxial and adaxial |
| phymac-ab | phymac-ad | - | - | *Phytelephas* leaf, abaxial and adaxial |
| ravmad-ab | ravmad-ad | - | - | *Ravenalia* leaf, abaxial and adaxial |


See Appendix II of the Concept note (https://docs.google.com/document/d/1W4qnylcvcscRP1e4GUldb6nxxLsuR12HL7QZ-eJv3zg/edit?usp=sharing) for instrument and measurement settings and measurement protocol (minimum five measurements per material, take white reference between materials). 

#### The simplest starting filename
```
[TC]<material>_<IDX>.<ext>
```
Examples: `fab2_0000.sig`, `TCfab2_0000.sed`

The `TC` prefix on the material is optional and will be stripped internally.

#### Full harmonization filename
```
PIdataHarmonization2026_HC<herbariumCode>_TC<material>_kit<n>_<IDX>.<ext>
```
Example: `PIdataHarmonization2026_HCHUH_TCfab2_kit1_0000.sig`

This is the format files will be converted to.

**WARNING:** The full filename convention does not, by itself, distinguish measurements made with **different foreoptics or optical setup configurations on the same instrument**. These differences can only be interpreted reliably if the files have distinct `measurementIndex` values and are linked to the metadata spreadsheet. If different optical setups are used, consider adding `SN<sessionId>` or another short foreoptic/optical setup identifier to the full filename (e.g., `SVC-8deg`). This is not a problem if the only difference is the instrument model because of different file extensions.

---

## Installing the spectra-harmonization GitHub repo

Terminal instructions:

Navigate to the folder where you want to clone the repository. For example:

```bash
cd ~/Documents/GitHub
```

Clone the repository from GitHub:

```bash
git clone https://github.com/YOUR-ORG-OR-USER/spectra-harmonization.git
```

Open the repository in RStudio by double-clicking `spectra-harmonization.Rproj`. This sets the working directory to the repository root, which all scripts assume.

Install required packages, if needed:

```r
install.packages(c("shiny", "shinyFiles", "tidyverse", "DT", "plotly", "fs","spectrolab","ggrepel"))
```

---

## Running the Shiny QC and metadata app (`spectra-files-qc-app.R`)

The Shiny app supports interactive visual QC, file flagging, metadata capture, and standardized file export.

### Run the app

```r
shiny::runApp("spectra-files-qc-app.R")
```

### Workflow

1. **Select a directory** containing your raw spectral files.
2. **Review file checks** — the app reports unrecognised filenames, materials with fewer than five measurements, and any of the 19 expected materials not present in the file set.
3. **Review QC plots** — spectra are plotted per material. Click a spectrum line or use the checkboxes to flag individual files for removal.
4. **Fill in the metadata form** — instrument, optical setup, operator, and measurement settings.
5. **Export** — unflagged files are copied to `shiny_qc_outputs/good_files_full_filenames/` using the full harmonization filename format. A metadata CSV is written alongside.

#### What constitutes a 'bad' spectral measurement that should be flagged?
Measurements with **differences in the shape of the spectral profile** should receive scrutiny and probably be removed. Differences in magnitude of the curves (ca. 5-10%) with no discernable shape difference are OK. Curves that look like **visual outliers**, however, should probably be removed.

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

Raw spectral files and outputs should generally not be committed to the repository. Add data file and outputs folders to `.gitignore` if you have moved such things into the repo folder.
