library(qtl2)
library(tidyverse)

## Load Data
#path <- "/home/c-hoyts/qtl_data/modular/"
path <- "C:/Users/c-hoyts/Documents/GitHub/jax-ssp19/QTL_mapping/modular/"
load(paste0(path, "DO_mESC_paired_eQTL_forMapping.RData"))
load(paste0(path,"peaks_for_Stephanie.RData"))
load(paste0(path, "ESC_eQTL_effects.RData"))
load(paste0(path, "scored_qtls.RData"))

## function to plot any peak given chromosome number and gene id
plot_CC <- function(chrom, gene, mgi, annotation = NULL){
  out <- scan1(probs, exprZ[, gene, drop=FALSE],
                   kinship_loco, addcovar=covar)
  eff <- scan1blup(probs[, chrom], exprZ[, gene, drop=FALSE],
                       kinship_loco[[chrom]], addcovar=covar)
  peak <- find_peaks(out, gmap, threshold=7) %>% filter(chr == chrom) %>% .$pos
  par(mar = c(8, 4, 4, 2) + 0.1)
  plot_coefCC(eff, map = gmap[chrom], scan1_output=out, main = mgi, legend="bottomleft")
  ## add vertical line at peak
  #abline(v = peak, col = "red", lty=2, lwd = 2)
  ## add info about scoring & clustering
  mtext(text = annotation, side = 1, outer = TRUE)
}

par(mfrow = c(2,2))
plot_CC(13, "ENSMUSG00000078994", "Zfp429") 
plot_CC(2, "ENSMUSG00000027597", "Ahcy")
plot_CC(12, "ENSMUSG00000021094", "Dhrs7")
plot_CC(17, "ENSMUSG00000052031", "Tagap1")





plot_CC(15, "ENSMUSG00000037280", "Galnt6")


###########################################################################################
## from Selcan

getplots <-function(g.name,probs,exprZ,kinship_loco,covar,p.chr, dist=FALSE,thr=7.0, prot=FALSE, gene.id, ...){
  out <- scan1(probs, exprZ[,gene.id , drop=FALSE],
               kinship_loco, addcovar=covar )
  peak <- find_peaks(out, gmap, drop=1.5,
                     threshold=thr, thresholdX=thr)
  if(p.chr %in% c(1:19) & dim(peak) > 0){
    peak <- peak %>% filter(chr == p.chr)
  }
  for(i in 1:dim(peak)[1]){
    chrom <- peak$chr[i]
    eff <- scan1blup(probs[,chrom], exprZ[,gene.id , drop=FALSE],
                     kinship_loco[[chrom]], addcovar=covar )
    # Let's find the genomic position - marker that is closest to my gene to abline
    peak.pos <- peak$pos[i]
   #  par(mfrow=c(2,1))
   #  plot(out, map = gmap, chr = chrom , gap=0, col="dark green", xlab="cM", bgcolor="gray95",
   #   main=paste0(g.name," (",gene.id,") ","plot"))
   # abline(h=thr,lwd=1)
   # abline(v=peak.pos,col="red",lwd=2,lty=3)
   if(peak.pos >30){plot_coefCC(eff, map = gmap[chrom], main=paste0(g.name," (",gene.id,") ","plot"),bgcolor="gray95",legend="bottomleft")}
   if(peak.pos <30){plot_coefCC(eff, map = gmap[chrom], main=paste0(g.name," (",gene.id,") ","plot"),bgcolor="gray95",legend="bottomright")}
   abline(v=peak.pos,col="red",lwd=2,lty=3)
 }
}

getplots("Zfp429", probs = probs, exprZ = exprZ, kinship_loco = kinship_loco, covar = covar, p.chr = 13, gene.id = "ENSMUSG00000078994")
