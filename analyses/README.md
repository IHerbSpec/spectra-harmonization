# Cross-instrument spectral characterization

Reproducible R analyses for the IHerbSpec harmonization project. The first question is whether spectra differ in **magnitude, shape, or both** across datasets, instruments, and optical configurations.

This workflow prepares spectra on a **400–2400 nm grid at 5 nm intervals**, plots means with central 90% percentile bands, and compares PCA before and after vector normalization. It also generates descriptive comparisons of confirmed shared physical targets. Constant-offset tests and transfer-function models are deferred.

The organization follows [herbarium-spectra](https://github.com/Erythroxylum/herbarium-spectra), with numbered analysis scripts, shared helpers, and separate plotting code. All new analysis code is under `analyses/`. A repository-level GitHub Actions workflow runs the numerical tests and synthetic example.

## Quick start

Clone the repository and open `spectra-harmonization.Rproj` in RStudio:

```bash
git clone https://github.com/IHerbSpec/spectra-harmonization.git
cd spectra-harmonization
```

From the **repository root** in R:

```r
# Run an end-to-end example using synthetic data.
source("analyses/R/00_make_example.R")
source("analyses/run_all.R")
```

Or in a terminal:

```bash
Rscript analyses/R/00_make_example.R
Rscript analyses/run_all.R
```

The example uses **base R only** and writes to `analyses/outputs/example/`. The 120 example spectra are synthetic, with deliberately introduced magnitude and shape differences. They are not project results.

An existing nonempty output directory causes an error. Choose a new `output_dir` or deliberately set `overwrite = TRUE` in your config to replace generated results. Check `RUN_STATUS.txt`: only `COMPLETE` indicates the full workflow finished.

## Repository structure

```text
analyses/
  README.md
  run_all.R                         # Run the complete workflow
  config/
    analysis.R                      # Range, spacing, FWHM, output settings
    datasets.example.csv            # Synthetic example manifest
    datasets.template.csv           # Starting point for your datasets
  R/
    00_make_example.R               # Generate synthetic example inputs
    01_import_preprocess.R          # Import, validate, Gaussian-resample
    02_characterize.R               # Summaries, PCA, matched comparisons
    03_plot.R                       # Write figures and run summary
  auxiliary/
    io.R                            # CSV and native-file adapters
    preprocessing.R                 # Gaussian resampling, normalization
    characterization.R              # Summary and comparison functions
  plotting/
    plotting_functions_spectra.R     # PDF and PNG figure functions
  data/
    raw/                            # Local native files or spectral CSVs
    metadata/template.csv           # Canonical metadata example
    example/                        # Generated synthetic inputs
  docs/
    methods.md                      # Resampling and interpretation
    input_formats.md                # Input schemas and QC-app mapping
  tests/
    run_tests.R                     # Numerical and import-contract checks
    check_example.R                 # End-to-end output checks
  outputs/                          # Generated locally, ignored by Git
```

## Use your own datasets

1. Place local inputs under `analyses/data/`, or use absolute paths in the manifest. Unzip any QC-app SHARE archives first.
2. Copy `config/datasets.template.csv` to `config/datasets.local.csv`. Add **one row per institution/instrument/optic/session dataset**. Split datasets when the wavelength grid differs.
3. Supply a per-spectrum metadata CSV. See [input formats](docs/input_formats.md). The QC app's `fullFilename`, `targetClass`, and `kitNumber` columns are accepted directly.
4. Copy `config/analysis.R` to `config/local.R`. Set `manifest = "config/datasets.local.csv"` and choose a new `output_dir`, such as `"outputs/2026-11-21"`.
5. For native files, install the optional reader: `install.packages("spectrolab")`. CSV inputs require no additional packages.
6. Run from the repository root:

```r
options(iherbspec.config = "config/local.R")
source("analyses/run_all.R")
```

Or:

```bash
Rscript analyses/run_all.R analyses/config/local.R
```

All paths **inside the config and manifest** are relative to `analyses/`, unless absolute. `run_all.R` also works when the working directory is `analyses/`. Numbered scripts share context through the runner; run the runner rather than sourcing stages independently.

No observations are excluded automatically for their appearance or reflectance magnitude. Address visual QC with the existing app first. Included data must have complete spectra, unique IDs, exact metadata matches, adequate wavelength coverage, and no source-band gaps above the configured limit. Failures stop the run with the dataset ID. Values outside [0,1] remain in the data and are reported for review.

## Analysis workflow

### 1. Import and preprocess

- Read each setup independently and join spectra to metadata by explicit ID.
- Convert percentages to fractions only when the manifest specifies `percent`.
- Retain native wavelength support for complete Gaussian windows.
- Resample onto `seq(400, 2400, by = 5)` using a Gaussian target response with **FWHM = 5 nm** by default.
- Export the resulting **401 bands**, metadata, QC table, and normalized spectra.

This follows the **target-FWHM Gaussian weighting** used by `CWT::resampling_FWHM()` in the earlier herbarium-spectra workflow. It is not linear interpolation. The output spacing and target Gaussian FWHM are separate settings. This operation does not guarantee equal effective resolution across instruments. See [methods](docs/methods.md).

### 2. Characterize magnitude and shape

- Mean and 5th–95th percentiles by dataset, kit, material, and physical sample when supplied.
- Mean reflectance and vector length for each spectrum.
- PCA of painted standards before and after vector normalization, using the same spectra and wavelength range.
- Matched comparisons across instrument families and across optical configurations within the **same physical SVC instrument**.
- Difference spectra, signed mean differences, RMSE, MAE, and spectral angle between dataset means on a shared physical target. These are descriptive summaries, not statistical tests.

When `physical_sample_id` is unavailable, overall summaries and PCA still run. Matched comparisons skip those observations and explain the omission. **Assign the same physical ID only to measurements of the same actual piece of material and surface**, not simply the same material type in different kits.

### 3. Plot

Figures use consistent dataset colors and include replicate counts. The same plotting code writes multipage PDFs and PNGs. Empty comparisons and insufficient PCA data generate `*_SKIPPED.txt` files.

| Figure prefix | Presentation use |
| --- | --- |
| `04_pnt_*` | Painted standards, raw reflectance and normalized values |
| `05_tvk_*` | Tyvek, raw reflectance and normalized values |
| `06_pca_reflectance` | Painted-standard PCA before normalization |
| `07_matched_instruments` | Matched physical targets across instrument families |
| `08_matched_svc_optics` | Matched targets across optics on the same SVC instrument |
| `10_pca_normalized` | Painted-standard PCA after normalization |
| `paired_mean_differences` | Wavelength-resolved differences between matched dataset means |

## Outputs

Each output directory contains:

| Folder | Contents |
| --- | --- |
| `processed/` | Reflectance and normalized CSVs, metadata, combined RDS, and fitted PCA objects |
| `tables/` | Mean/percentile spectra, inventory, QC flags, magnitude summaries, PCA results, matched sample lists, comparison metrics, and dataset colors |
| `figures/` | Presentation figures in PDF and PNG formats |
| `provenance/` | Config, dataset manifest, input/code MD5 checksums, method settings, and R/package session information |

Raw files remain unchanged. Local observations, local configs, and outputs are ignored by Git. Version your final config and manifest with the archived dataset when sharing a reproducible release. Native-file reader versions are recorded in `sessionInfo.txt`.

## Interpretation

- Percentile bands describe variation among scans. They are **not confidence intervals**.
- Kit/material/physical-sample groups stay separate in mean plots. Multiple scans do not imply multiple independent physical samples.
- PCA centers wavelength columns and does not scale them to unit variance. Axes from separately fitted raw and normalized PCAs need not align.
- Instrument comparisons retain dataset labels and optical configurations. Instrument, institution, operator, and setup effects may remain confounded.
- Normalization removes vector magnitude. Remaining differences motivate subsequent harmonization work.
- A smaller discrepancy is useful only if meaningful differences among materials remain. Formal validation and transfer models are later steps.

## Tests

```bash
Rscript analyses/tests/run_tests.R
Rscript analyses/R/00_make_example.R
Rscript analyses/run_all.R
Rscript analyses/tests/check_example.R
```

The numerical tests check Gaussian convolution against an analytical result, constant preservation, endpoint support, gaps, normalization, metadata order, units, matching rules, input errors, and PCA. If `CWT` is installed, they also compare Gaussian output directly against `CWT::resampling_FWHM()`. The GitHub Actions workflow runs the base-R checks and the full synthetic example.

**Initial delivery status:** the R tests and example have not run in the authoring environment. See [validation status](docs/validation_status.md) before using this first version.
