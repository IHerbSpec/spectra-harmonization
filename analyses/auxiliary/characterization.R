# Exploratory summaries only. Technical scans are not independent materials.
group_indices <- function(m, cols) {
  split(seq_len(nrow(m)),do.call(paste,c(m[cols],sep="\r")))
}

spectral_summary <- function(x,m,wavelengths) {
  cols <- c("dataset_id","kit_id","material_id","material_group","physical_sample_id","instrument_family",
            "instrument_model","instrument_id","optic_id")
  groups <- group_indices(m,cols)
  do.call(rbind,lapply(groups,function(ii) {
    y <- x[ii,,drop=FALSE]
    q <- apply(y,2,quantile,probs=c(0.05,0.95),names=FALSE,type=7)
    cbind(m[rep(ii[1],length(wavelengths)),cols,drop=FALSE],
          data.frame(wavelength_nm=wavelengths,n_spectra=length(ii),mean=colMeans(y),
                     q05=q[1,],q95=q[2,],sd=if(length(ii)>1L) apply(y,2,sd) else NA_real_))
  }))
}

fit_pca <- function(x,m) {
  if (nrow(x)<3L) return(list(skipped="PCA needs at least 3 spectra"))
  keep <- apply(x,2,var)>1e-20
  if (sum(keep)<2L) return(list(skipped="PCA needs at least 2 variable bands"))
  fit <- prcomp(x[,keep,drop=FALSE],center=TRUE,scale.=FALSE,rank.=min(10L,nrow(x)-1L,sum(keep)))
  if (sum(fit$sdev>max(fit$sdev)*1e-10)<2L) return(list(skipped="PCA matrix has rank below 2"))
  list(model=fit,scores=cbind(m,as.data.frame(fit$x)),
       loadings=data.frame(wavelength_nm=as.numeric(colnames(x)[keep]),fit$rotation,check.names=FALSE),
       variance=data.frame(component=paste0("PC",seq_along(fit$sdev)),variance_fraction=fit$sdev^2/sum(fit$sdev^2)))
}

pairwise_means <- function(x,m,wavelengths) {
  eligible <- !is.na(m$physical_sample_id)
  x <- x[eligible,,drop=FALSE];m <- m[eligible,,drop=FALSE]
  empty <- list(metrics=data.frame(),differences=data.frame())
  if (!nrow(m)) return(empty)
  groups <- group_indices(m,c("kit_id","material_id","physical_sample_id"))
  metrics <- differences <- list(); counter <- 0L
  for (ii in groups) {
    setups <- split(ii,m$dataset_id[ii])
    if (length(setups)<2L) next
    means <- lapply(setups,function(j) colMeans(x[j,,drop=FALSE]))
    pairs <- combn(names(setups),2,simplify=FALSE)
    for (p in pairs) {
      a <- m[setups[[p[1]]][1],,drop=FALSE];b <- m[setups[[p[2]]][1],,drop=FALSE]
      av <- means[[p[1]]];bv <- means[[p[2]]];delta <- av-bv
      denom <- sqrt(sum(av^2)*sum(bv^2))
      angle <- if (denom<=0) NA_real_ else acos(max(-1,min(1,sum(av*bv)/denom)))*180/pi
      kind <- if (a$instrument_id!=b$instrument_id) "between_instruments" else
        if (a$optic_id!=b$optic_id) "within_instrument_optics" else "between_datasets_same_setup"
      info <- data.frame(kit_id=a$kit_id,material_id=a$material_id,material_group=a$material_group,
        physical_sample_id=a$physical_sample_id,dataset_a=p[1],dataset_b=p[2],
        family_a=a$instrument_family,family_b=b$instrument_family,
        instrument_a=a$instrument_id,instrument_b=b$instrument_id,
        optic_a=a$optic_id,optic_b=b$optic_id,comparison=kind,
        n_a=length(setups[[p[1]]]),n_b=length(setups[[p[2]]]))
      counter <- counter+1L
      metrics[[counter]] <- cbind(info,data.frame(mean_a=mean(av),mean_b=mean(bv),
        signed_mean_difference=mean(delta),rmse=sqrt(mean(delta^2)),mae=mean(abs(delta)),spectral_angle_deg=angle))
      differences[[counter]] <- cbind(info[rep(1,length(wavelengths)),,drop=FALSE],
        data.frame(wavelength_nm=wavelengths,difference_a_minus_b=delta))
    }
  }
  if (!counter) return(empty)
  list(metrics=do.call(rbind,metrics),differences=do.call(rbind,differences))
}

matched_rows <- function(m,mode) {
  ok <- !is.na(m$physical_sample_id)
  if (mode=="svc_optics") ok <- ok & toupper(m$instrument_family)=="SVC"
  ii <- which(ok)
  if (!length(ii)) return(integer())
  cols <- c("kit_id","material_id","physical_sample_id")
  if (mode=="svc_optics") cols <- c(cols,"instrument_id")
  groups <- group_indices(m[ii,,drop=FALSE],cols)
  unlist(lapply(groups,function(j) {
    actual <- ii[j]
    n <- if (mode=="svc_optics") length(unique(m$optic_id[actual])) else length(unique(m$instrument_family[actual]))
    if (n>=2L) actual else integer()
  }),use.names=FALSE)
}
