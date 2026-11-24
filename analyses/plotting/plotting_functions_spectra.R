# Base-R figures keep the CSV/example workflow free of package dependencies.
# PDF pages and 180-dpi PNGs use the same drawing functions.
safe_name <- function(x) gsub("[^A-Za-z0-9_-]","_",x)
draw_to_file <- function(path,draw,width=12,height=7) {
  if (grepl("[.]pdf$",path)) pdf(path,width=width,height=height,onefile=TRUE) else
    png(path,width=width,height=height,units="in",res=180)
  on.exit(dev.off(),add=TRUE)
  draw()
}
plot_note <- function(path,message) {
  writeLines(message,paste0(path,"_SKIPPED.txt"))
}
legend_panel <- function(labels,colors) {
  par(mar=c(0,0,0,0));plot.new()
  legend("center",legend=labels,col=colors,lty=1,lwd=2,bty="n",ncol=min(3,length(labels)),cex=0.8)
}

plot_mean_pages <- function(summary,path,colors,y_label="Reflectance (fraction)",title_prefix="") {
  if (!nrow(summary)) return(plot_note(path,"No matching spectra"))
  panels <- group_indices(summary,c("kit_id","material_id","physical_sample_id"))
  limits <- range(summary$q05,summary$q95,finite=TRUE)
  if (diff(limits)==0) limits <- limits+c(-0.01,0.01)
  draw_panel <- function(ii) {
    d <- summary[ii,,drop=FALSE]
    layout(matrix(1:2,2,1),heights=c(5,1.1))
    par(mar=c(4.3,4.8,3.5,1))
    physical <- if(is.na(d$physical_sample_id[1])) "physical sample unspecified" else d$physical_sample_id[1]
    plot(range(d$wavelength_nm),limits,type="n",xlab="Wavelength (nm)",ylab=y_label,
         main=paste0(title_prefix,d$material_id[1]," / kit ",d$kit_id[1]," / ",physical))
    groups <- split(seq_len(nrow(d)),d$dataset_id)
    for (id in names(groups)) {
      dd <- d[groups[[id]],,drop=FALSE];dd <- dd[order(dd$wavelength_nm),]
      polygon(c(dd$wavelength_nm,rev(dd$wavelength_nm)),c(dd$q05,rev(dd$q95)),
              col=adjustcolor(colors[[id]],alpha.f=0.16),border=NA)
      lines(dd$wavelength_nm,dd$mean,col=colors[[id]],lwd=1.7)
    }
    labels <- vapply(names(groups),function(id) paste0(id," (n=",d$n_spectra[groups[[id]][1]],")"),character(1))
    legend_panel(labels,colors[names(groups)])
  }
  draw_to_file(paste0(path,".pdf"),function() for(ii in panels) draw_panel(ii))
  for (k in seq_along(panels)) draw_to_file(paste0(path,"_",sprintf("%02d",k),".png"),function() draw_panel(panels[[k]]))
}

plot_pca <- function(fit,path,colors,title) {
  if (!is.null(fit$skipped)) return(plot_note(path,fit$skipped))
  draw <- function() {
    d <- fit$scores
    materials <- sort(unique(d$material_id))
    shapes <- setNames(rep(c(16,17,15,18,1,2,0,5),length.out=length(materials)),materials)
    layout(matrix(c(1,2,3,3),2,2,byrow=TRUE),heights=c(5,1.2))
    par(mar=c(4.5,4.5,3,1))
    vari <- fit$variance$variance_fraction
    plot(d$PC1,d$PC2,pch=shapes[d$material_id],col=colors[d$dataset_id],
         xlab=sprintf("PC1 (%.1f%%)",100*vari[1]),ylab=sprintf("PC2 (%.1f%%)",100*vari[2]),main=title)
    legend("topright",legend=materials,pch=shapes,col="grey30",bty="n",cex=0.7)
    matplot(fit$loadings$wavelength_nm,fit$loadings[c("PC1","PC2")],type="l",lty=1:2,col=c("#1B6A54","#995A35"),
            xlab="Wavelength (nm)",ylab="Loading",main="PCA loadings")
    abline(h=0,col="grey75");legend("topright",c("PC1","PC2"),lty=1:2,col=c("#1B6A54","#995A35"),bty="n")
    ids <- sort(unique(d$dataset_id));legend_panel(ids,colors[ids])
  }
  draw_to_file(paste0(path,".pdf"),draw)
  draw_to_file(paste0(path,".png"),draw)
}

plot_differences <- function(d,path) {
  if (!nrow(d)) return(plot_note(path,"No explicitly matched physical samples across datasets"))
  groups <- group_indices(d,c("kit_id","material_id","physical_sample_id","dataset_a","dataset_b"))
  draw_panel <- function(ii) {
    z <- d[ii,,drop=FALSE];z <- z[order(z$wavelength_nm),]
    par(mar=c(5,5,5,2))
    plot(z$wavelength_nm,z$difference_a_minus_b,type="l",lwd=1.5,col="#1B6A54",
         xlab="Wavelength (nm)",ylab="Difference in mean reflectance (A - B)",
         main=paste0(z$material_id[1]," / ",z$physical_sample_id[1],"\n",z$dataset_a[1]," minus ",z$dataset_b[1]))
    abline(h=0,lty=2,col="grey60")
  }
  draw_to_file(paste0(path,".pdf"),function() for(ii in groups) draw_panel(ii))
  for (k in seq_along(groups)) draw_to_file(paste0(path,"_",sprintf("%03d",k),".png"),function() draw_panel(groups[[k]]))
}
