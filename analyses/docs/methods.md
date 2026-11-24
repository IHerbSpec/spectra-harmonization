# Preprocessing and characterization methods

## Range, spacing, and FWHM

The default output is `seq(400, 2400, 5)`: 401 bands, including both endpoints.

The Gaussian response at output wavelength c is proportional to:

```text
sigma = FWHM / (2 * sqrt(2 * log(2)))
weight(lambda, c) = exp(-0.5 * ((lambda - c) / sigma)^2)
R_out(c) = sum(weight(lambda, c) * R_in(lambda)) / sum(weight(lambda, c))
```

This independently implemented discrete channel weighting follows the algorithm in [CWT's resampling_FWHM.cpp](https://github.com/Antguz/CWT/blob/29d626a126f244ddb3a5cc6faf18c376723adcf4/src/resampling_FWHM.cpp). The earlier [herbarium-spectra/R/spectra_transform.R](https://github.com/Erythroxylum/herbarium-spectra/blob/20ed89324307df404d1bc70050d21a22f6f17f78/R/spectra_transform.R) calls that function with a 5 nm output interval and a constant target FWHM of 5 nm. The new workflow uses explicit wavelength values, rather than positional column slices.

Both parameters default to 5 nm, but they mean different things:

- **Spacing** locates the output bands.
- **Target Gaussian FWHM** sets the resampling kernel width.
- **Native instrument FWHM** describes the measurement's existing spectral response and is not inferred from spacing.

This is Gaussian target-band resampling, **not native-resolution matching**. Applying the same kernel to measurements from different native spectral responses does not make those responses identical or recover lost resolution. Native-FWHM matching would require reliable response information and a separate, explicit broadening plan. Do not describe this initial preprocessing as resolving instrument differences.

The discrete weights also follow the source channel density. On strongly irregular grids they are not equivalent to wavelength-width-weighted numerical integration. Record this choice and examine unusual native grids before combining them. Maximum source gaps are checked over the working range and kernel margins. A denser grid created upstream cannot establish finer native resolution.

## Endpoint handling

For each output wavelength, native coverage must extend at least 3 Gaussian standard deviations on each side. For FWHM = 5 nm this is approximately 6.37 nm, so the endpoints require measurements to approximately 393.63 and 2406.37 nm.

The Gaussian uses all available input channels, as in CWT, after checking that coverage condition. The tiny Gaussian tails outside source coverage are consequently truncated and weights renormalized. Retain native bands outside 400–2400 for this calculation. The **final output** alone is restricted to 400–2400. Do not pre-trim a broad native spectrum to exactly those endpoints.

Datasets lacking coverage fail rather than extrapolate, zero-fill, or silently change the comparison range. Exclude a short-range dataset explicitly with `include=FALSE` and an `exclusion_reason`, or run a separate, documented common-overlap configuration. Native gaps above `max_source_gap_nm` also fail.

## Vector normalization

For each resampled spectrum x, divide by `sqrt(sum(x^2))` across the full configured wavelength interval. Save that norm. Normalized values are dimensionless and should not be labeled reflectance. Reject zero vectors. No centering, standard-normal-variate transformation, additive-offset correction, or transfer model is applied to individual spectra.

## Means and 90% bands

Group by dataset, kit, material, physical sample (where supplied), and instrument configuration. Calculate arithmetic means, sample SD, and the 5th and 95th percentiles using R's type-7 quantiles. Report n scans. With only five scans, these percentile estimates are descriptive and coarse. Never label them confidence intervals. When physical IDs are missing, the kit/material group may combine pieces of the same material; review the input metadata before interpreting it as technical variation.

## PCA

Use all individual painted-standard scans, center each band, and do not scale bands to unit variance. Fit raw-reflectance and normalized PCAs separately on the same observations and band range. Drop bands with variance <= 1e-20, record retained bands in loadings, and skip matrices with insufficient rank. Retain up to 10 PCs, while variance fractions use total spectral variance. Compare grouping patterns rather than treating axes of separate PCAs as directly interchangeable.

## Matched comparisons

`physical_sample_id`, `kit_id`, and `material_id` identify eligible shared targets. Separate surfaces must have separate material/physical IDs. Average scans within each dataset and compare those means, without pairing arbitrary replicate indices.

For A minus B, export wavelength-specific differences, their signed mean, RMSE, MAE, and spectral angle in degrees. Dataset IDs determine A/B order, and the tables record that direction. Differences are descriptive. No hypothesis tests, confidence intervals, constant-offset fits, or transfer-function validation appear in this phase.

The instrument figure requires at least two instrument families on the same target. The SVC optic figure requires at least two optic IDs on the same target **and physical instrument ID**. Paired tables additionally retain between-instrument comparisons within the same family and between-dataset comparisons of the same setup. No institution or instrument is assumed to be the true reflectance reference.

## Later extensions

After inspecting magnitude and shape differences, specify evaluation targets and independent validation units before fitting transformations. Potential next methods include native-resolution matching, EPO, PLSR transfer, and an explicit strategy for restoring reflectance magnitude. Constant-offset testing remains deferred.
