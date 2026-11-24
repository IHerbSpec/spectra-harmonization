# Base-R numerical and input-contract tests. Run from repository root.
local({
  root <- if(file.exists("analyses/auxiliary/io.R")) "analyses" else "."
  for (f in c("io.R","preprocessing.R","characterization.R")) source(file.path(root,"auxiliary",f),local=TRUE)
  n_tests <- 0L
  check <- function(ok,label) {if(!isTRUE(ok)) stop("FAIL: ",label);n_tests <<- n_tests+1L}
  fails <- function(code) inherits(tryCatch({force(code);NULL},error=identity),"error")
  w <- seq(390,2410,by=1);target <- seq(400,2400,5)
  const <- matrix(0.4,2,length(w))
  y <- gaussian_resample(const,w,target)
  check(all(dim(y)==c(2,401)),"401 output bands")
  check(max(abs(y-0.4))<1e-12,"Gaussian weights preserve a constant")
  slope <- matrix(w/2500,1)
  z <- gaussian_resample(slope,w,target)
  check(max(abs(z-as.numeric(target)/2500))<1e-6,"Linear spectrum and endpoint support")
  # Analytical convolution of a Gaussian feature with a Gaussian kernel.
  width <- 12;sigma <- 5/(2*sqrt(2*log(2)))
  signal <- matrix(exp(-0.5*((w-1500)/width)^2),1)
  actual <- gaussian_resample(signal,w,target)
  expected <- width/sqrt(width^2+sigma^2)*exp(-0.5*((target-1500)/sqrt(width^2+sigma^2))^2)
  check(max(abs(actual-expected))<1e-10,"Gaussian feature convolution matches analytic result")
  check(fails(gaussian_resample(const[,11:2011],w[11:2011],target)),"Reject premature trimming")
  check(fails(gaussian_resample(const[,c(1:500,520:length(w))],w[c(1:500,520:length(w))],target)),"Reject spectral gaps")
  check(fails(gaussian_resample(const[,c(1,1,2:length(w))],w[c(1,1,2:length(w))],target)),"Reject duplicate wavelengths")
  v <- vector_normalize(rbind(c(1,2,3),c(4,8,12)))
  check(max(abs(v$spectra[1,]-v$spectra[2,]))<1e-12,"Normalization removes multiplicative scale")
  check(max(abs(sqrt(rowSums(v$spectra^2))-1))<1e-12,"Normalized vector length is one")
  check(fails(vector_normalize(matrix(0,1,3))),"Reject zero spectrum")
  m <- data.frame(dataset_id=c("a","a","b","b"),kit_id="1",material_id="pnt1",material_group="pnt",
    physical_sample_id="kit1_pnt1",instrument_family=c("SVC","SVC","ASD","ASD"),
    instrument_model=c("s","s","a","a"),instrument_id=c("s1","s1","a1","a1"),optic_id="contact")
  x <- rbind(c(1,2,3),c(3,4,5),c(2,3,4),c(4,5,6));colnames(x)<-c("400","405","410")
  sm <- spectral_summary(x,m,c(400,405,410))
  check(all(sm$n_spectra==2),"Summary replicate counts")
  check(abs(sm$q05[1]-1.1)<1e-12 && abs(sm$q95[1]-2.9)<1e-12,"Type-7 90 percent interval")
  pairs <- pairwise_means(x,m,c(400,405,410))
  check(nrow(pairs$metrics)==1 && abs(pairs$metrics$rmse-1)<1e-12,"Pairwise mean RMSE")
  check(abs(pairs$metrics$signed_mean_difference+1)<1e-12,"Difference direction A minus B")
  m$physical_sample_id[3:4] <- "different_target"
  check(nrow(pairwise_means(x,m,c(400,405,410))$metrics)==0,"Never match different physical targets")
  check(length(matched_rows(m,"instruments"))==0,"No instrument match on different targets")
  m$physical_sample_id <- "shared";m$instrument_family <- "SVC";m$instrument_id <- "s1";m$optic_id <- c("a","a","b","b")
  check(length(matched_rows(m,"svc_optics"))==4,"Within-instrument optic match")
  m$instrument_id[3:4] <- "s2"
  check(length(matched_rows(m,"svc_optics"))==0,"Optic comparison requires same instrument")
  set.seed(19);p <- fit_pca(matrix(rnorm(10*20),10,dimnames=list(NULL,as.character(1:20))),data.frame(sample_key=1:10))
  check(is.null(p$skipped) && abs(sum(p$variance$variance_fraction)-1)<1e-12,"PCA variance fractions")
  check(!is.null(fit_pca(matrix(1,2,10),data.frame(id=1:2))$skipped),"Small PCA skips explicitly")
  # Wide/long import, percent conversion, metadata joins, and ID failures.
  tmp <- tempfile("iherbspec-tests-");dir.create(tmp);on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
  metadata <- data.frame(spectrum_id=c("b","a"),material_id=c("tvk","pnt1"),kit_id="1",physical_sample_id=c("kt","kp"))
  write_csv(metadata,file.path(tmp,"metadata.csv"))
  wide <- data.frame(spectrum_id=c("a","b"),rbind(rep(20,length(w)),rep(80,length(w))),check.names=FALSE)
  names(wide)[-1] <- as.character(w);write_csv(wide,file.path(tmp,"wide.csv"))
  row <- data.frame(dataset_id="test",format="wide_csv",spectra_path="wide.csv",metadata_path="metadata.csv",institution="test",
    instrument_family="SVC",instrument_model="test",instrument_id="s1",optic_id="contact",reflectance_units="percent")
  cfg <- list(wavelength_min_nm=400,wavelength_max_nm=2400,spacing_nm=5,target_fwhm_nm=5,max_source_gap_nm=5)
  loaded <- load_dataset(row,tmp,cfg)
  check(identical(loaded$metadata$material_id,c("pnt1","tvk")),"Join by ID rather than metadata row order")
  check(max(abs(loaded$spectra[1,]-0.2))<1e-12,"Explicit percent to fraction conversion")
  long <- data.frame(spectrum_id=rep(c("a","b"),each=length(w)),wavelength_nm=rep(w,2),reflectance=c(rep(20,length(w)),rep(80,length(w))))
  write_csv(long,file.path(tmp,"long.csv"));row$format<-"long_csv";row$spectra_path<-"long.csv"
  check(max(abs(load_dataset(row,tmp,cfg)$spectra-loaded$spectra))<1e-12,"Wide and long inputs agree")
  write_csv(rbind(long,long[1,]),file.path(tmp,"long.csv"))
  check(fails(read_long(file.path(tmp,"long.csv"))),"Reject duplicate long sample/wavelength rows")
  long <- long[-1,];write_csv(long,file.path(tmp,"long.csv"))
  check(fails(read_long(file.path(tmp,"long.csv"))),"Reject incomplete spectral grids")
  metadata$spectrum_id[1]<-"unmatched";write_csv(metadata,file.path(tmp,"metadata.csv"));row$format<-"wide_csv";row$spectra_path<-"wide.csv"
  check(fails(load_dataset(row,tmp,cfg)),"Reject unmatched metadata IDs")
  if(requireNamespace("CWT",quietly=TRUE)) {
    yy <- as.matrix(CWT::resampling_FWHM(spectra=signal,wavelengths=w,new_wavelengths=target,FWHM=rep(5,length(target))))
    check(max(abs(yy-actual))<1e-12,"Numerical agreement with CWT::resampling_FWHM")
  }
  # Optional native-reader smoke check on an existing repository .sed file.
  native_dir <- file.path(root,"..","test_data")
  if(requireNamespace("spectrolab",quietly=TRUE) && file.exists(file.path(native_dir,"pnt1_00099.sed"))) {
    nm <- data.frame(spectrum_id="pnt1_00099.sed")
    nd <- read_native(native_dir,nm)
    check(nrow(nd$spectra)==1L && length(nd$wavelengths)==ncol(nd$spectra),"Native file retains bands and sample identity")
    check(all(is.finite(gaussian_resample(nd$spectra,nd$wavelengths,target))),"Native .sed Gaussian preprocessing")
  }
  message("PASS: ",n_tests," checks")
})
