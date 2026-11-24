# Validation status of this initial delivery

The implementation includes numerical tests and a complete synthetic example, but **R execution has not been verified in the authoring environment**. R/Rscript was unavailable, and the environment could not install an R runtime. The GitHub connection rejected writes, so the supplied Actions workflow has not run either.

Before interpreting any outputs, run from the repository root:

```bash
Rscript analyses/tests/run_tests.R
Rscript analyses/R/00_make_example.R
Rscript analyses/run_all.R
Rscript analyses/tests/check_example.R
```

CSV inputs and the synthetic example use base R. Install `spectrolab` to exercise the optional native-file smoke test on the existing repository test data. Install `CWT` to exercise the optional direct comparison with `resampling_FWHM`.

The example and tests contain no project observations. Native ASD/SVC/SE data should be checked against their source software for units and reader behavior before scientific use.
