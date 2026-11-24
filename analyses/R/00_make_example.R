# Deterministic synthetic inputs for demonstrating the complete pipeline.
# No project observations are included or implied.
local({
  root <- if(file.exists("analyses/run_all.R")) "analyses" else if(file.exists("run_all.R")) "." else stop("Run from the repository root or analyses/")
  path <- file.path(root,"data/example")
  dir.create(path,recursive=TRUE,showWarnings=FALSE)
  set.seed(20261121)
  materials <- c("pnt1","pnt2","pnt3","pnt4","tvk","fab2")
  datasets <- c("demo_naturaspec","demo_svc_contact","demo_svc_lens","demo_asd")
  for (k in seq_along(datasets)) {
    w <- seq(390,2410,by=if(k==1) 2 else 1)
    material <- rep(materials,each=5)
    ids <- sprintf("synthetic_%03d",seq_along(material))
    x <- t(vapply(seq_along(ids),function(i) {
      j <- match(material[i],materials)
      level <- c(0.12,0.28,0.46,0.65,0.82,0.4)[j]
      curve <- level+0.025*sin(w/300)+0.015*(w-1400)/1000
      if (j>=5) curve <- curve-0.08*exp(-((w-1720)/100)^2)-0.13*exp(-((w-2250)/80)^2)
      scale <- c(1,0.97,1.04,1.02)[k]
      curve*scale+(k-1)*0.004*sin(w/180)+rnorm(1,0,0.002)+rnorm(length(w),0,0.0008)
    },numeric(length(w))))
    if(k==4) x <- x*100
    colnames(x) <- as.character(w)
    write.csv(data.frame(spectrum_id=ids,x,check.names=FALSE),file.path(path,paste0(datasets[k],"_spectra.csv")),row.names=FALSE)
    m <- data.frame(spectrum_id=ids,material_id=material,kit_id="DEMO1",
      material_group=ifelse(grepl("^pnt",material),"pnt",material),
      physical_sample_id=paste0("DEMO1_",material),replicate=rep(1:5,length(materials)))
    write.csv(m,file.path(path,paste0(datasets[k],"_metadata.csv")),row.names=FALSE)
  }
  message("Synthetic example files written to ",path)
})
