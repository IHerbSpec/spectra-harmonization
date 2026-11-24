# Input formats

## Dataset manifest

One row represents one dataset with a single instrument/optic configuration and wavelength grid. Use a different dataset ID for different sessions or grids. Inputs are explicit files or directories, not discovered from institution names.

| Column | Meaning |
| --- | --- |
| `dataset_id` | Unique stable ID, using letters, digits, underscores, or hyphens |
| `format` | `wide_csv`, `long_csv`, or `native` |
| `spectra_path` | CSV path, or directory of native files |
| `metadata_path` | Per-spectrum metadata CSV |
| `institution` | Institution or project label |
| `instrument_family` | Prefer `SVC`, `NaturaSpec`, or `ASD`; other labels are allowed |
| `instrument_model` | Actual instrument model |
| `instrument_id` | Globally unique physical unit ID; use the same ID across its optics |
| `optic_id` | Stable optical-configuration ID, including geometry where relevant |
| `reflectance_units` | Explicit `fraction` or `percent`; never guessed |
| `include` | `TRUE` or `FALSE` |
| `exclusion_reason` | Required when `include=FALSE`, e.g. superseded measurements or insufficient range |

Use separate manifest rows for lens or probe configurations. Put details such as field of view, geometry, operator, serial number, reference standard, and measurement area in the per-spectrum metadata. Additional metadata columns survive import. Do not place manifest-reserved columns in metadata: dataset_id, institution, instrument_family, instrument_model, instrument_id, optic_id, reflectance_units, sample_key.

## Per-spectrum metadata

```csv
spectrum_id,material_id,material_group,kit_id,physical_sample_id,replicate
pnt1_0001.sig,pnt1,pnt,1,kit1_pnt1,1
pnt1_0002.sig,pnt1,pnt,1,kit1_pnt1,2
tvk_0001.sig,tvk,tvk,1,kit1_tvk,1
```

Required: `spectrum_id`, `material_id`, `kit_id`. IDs must match the spectra exactly and uniquely within each dataset. IDs are read as character strings to retain leading zeroes. `sample_key` in outputs combines dataset ID and spectrum ID, so filenames may recur across different datasets.

`material_group` defaults to `pnt` for IDs beginning `pnt`, and otherwise to the material ID. You may set it explicitly for new materials. The runner uses `config$material_groups` for mean plots and `config$pca_material_group` for PCA.

`physical_sample_id` is optional for exploratory means/PCA and required for matched comparisons. The same ID must mean the same physical object and surface across datasets. It must map to only one kit and material. Different kits containing the same paint are different physical targets. Leave this field blank when identity cannot be established.

## Reuse metadata from the Shiny QC app

The following aliases are automatic:

| QC export column | Analysis column |
| --- | --- |
| `fullFilename` | `spectrum_id` |
| `targetClass` | `material_id` |
| `kitNumber` | `kit_id` |

The other QC fields remain available, including `instrumentModel`, `opticalSetupDescription`, `lensFieldOfView`, `measurementAreaDiameter`, and `comment`. Supply the manifest's instrument/optic IDs explicitly. Do not infer an instrument identity or geometry from file extension alone.

For matched comparisons, add and verify `physical_sample_id`. Split the metadata by setup if an export combines configurations. Use the exported renamed good files so `fullFilename` matches the files on disk.

## Wide CSV

One row per scan, with `spectrum_id` plus numeric wavelength headers in **nm**:

```csv
spectrum_id,390,391,392,393
scan01,0.14,0.14,0.15,0.15
scan02,0.15,0.15,0.16,0.16
```

This shortened example illustrates column structure only. Actual data must span the full output range plus Gaussian margins. Headers such as `X400` are accepted as well as `400` and `400.5`. Metadata belongs in the separate metadata CSV. Do not include a saved row-number column.

## Long CSV

```csv
spectrum_id,wavelength_nm,reflectance
scan01,390,0.14
scan01,391,0.14
scan02,390,0.15
scan02,391,0.15
```

Each sample/wavelength pair must occur exactly once. All scans in a dataset must share the same complete wavelength grid. Missing reflectance values stop the analysis. To adapt a custom format, export this schema or the wide schema rather than modifying the analysis stages.

## Native files

Set `format=native`, point `spectra_path` to a directory, and use exact filenames with extensions for `spectrum_id`. Only files explicitly listed in metadata are read. Other files in the directory are ignored. Supported reader extensions are `.sig`, `.sed`, `.asd`, and `.txt`, through `spectrolab::read_spectra(type="target_reflectance")`.

Each file must contain one reflectance spectrum. Reading each file individually avoids dependence on directory order and ASD sample-label conventions. Files within one manifest row must share a native wavelength grid. The analysis does not splice detectors, fit reference corrections, or alter manufacturer settings. Native formats unsupported by your installed reader should be converted to a documented reflectance CSV first.

Verify the native reader's output units against a known file and set `reflectance_units` accordingly. Do not assume all reader/instrument combinations use the same scale.
