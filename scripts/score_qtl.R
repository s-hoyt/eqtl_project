library(qtl2)
library(tidyverse)

## Load Data
args = commandArgs(trailingOnly = TRUE)
for(i in 2:(length(args)-1)){
  load(paste0(args[1], args[i]))
}
#path <- "/home/c-hoyts/qtl_data/"
# path <- "C:/Users/c-hoyts/Documents/QTL_mapping/"
path <- "C:/Users/c-hoyts/Documents/GitHub/jax-ssp19/QTL_mapping/modular/"
load(paste0(path, "ranking_save_tads.RData"))


qtl_score <- cbind(qtl_score, peaks_local_over7 %>% dplyr::select(peak, genes))

#version 2 - distance subtracted, and with snp threshold
score2_func <- function(qtls, thr = 150, snp_w = 1000, lod_w = 100, biallelic_w = 10, distance_w = 10000, biotype_w = 10, tad_w = 10){
  num_snps <- qtls["num_snps"] %>% unlist() %>% as.integer()
  if(is.na(num_snps)){
    return(NA)
  }
  if(num_snps == 0 | num_snps == -1){
    snps <- -(thr * snp_w)
  }else if(num_snps > thr){
    snps <- -(num_snps/thr) * snp_w
  }else{
    snps <- (1/num_snps) * snp_w
  }
  tad <- 0
  if(!(is.na(qtls["peak"]) | is.na(qtls["genes"]))){
    if((qtls["peak"] %>% as.numeric) == (qtls["genes"] %>% as.numeric)){
      tad <- tad_w
    }
  }
  lod <- (qtls["lod"] %>% as.double()) * lod_w
  dist <- (qtls["peak2gene"] %>% as.double) / distance_w
  biallelic <- biallelic_w
  pseudo_biallelic <- (!(qtls["biallelic"] %>% as.logical())) & (qtls["num_snps"] != -1)
  if(!(qtls["biallelic"] %>% as.logical())){
    biallelic <- -biallelic_w
  }else if(pseudo_biallelic){
    biallelic <- .5*biallelic_w
  }
  if(qtls["gene_biotype"] == "protein_coding"){
    biotype <- biotype_w
  }else{
    biotype <- 0
  }
  score <- (lod + snps + biallelic + biotype + tad) - dist
  return(score)
}

score2 <- qtl_score %>% t() %>% as.data.frame() %>% purrr::map(score2_func, thr = 150, snp_w = 1000, lod_w = 100, biallelic_w = 10, 
                                                               distance_w = 10000, biotype_w = 10, tad_w = 10)
qtl_score <- qtl_score %>% mutate(score = score2)

## function to plot any peak given chromosome number and gene id, and optional annotation (ie scores, clusters)
plot_CC <- function(chrom, gene, annotation = NULL){
  out <- scan1(probs, exprZ[, gene, drop=FALSE],
               kinship_loco, addcovar=covar)
  eff <- scan1blup(probs[, chrom], exprZ[, gene, drop=FALSE],
                   kinship_loco[[chrom]], addcovar=covar)
  peak <- find_peaks(out, gmap, threshold=7) %>% filter(chr == chrom) %>% .$pos
 # par(mar = c(8, 4, 4, 2) + 0.1)
  plot_coefCC(eff, map = gmap[chrom], scan1_output=out, main = gene)
  ## add vertical line at peak
  abline(v = peak, col = "red", lty=2, lwd = 2)
  ## add info about scoring & clustering
  #mtext(text = annotation, side = 1, outer = TRUE)
}

#plot_CC(2, "ENSMUSG00000027014")

save(qtl_score, file = paste0(args[1], args[length(args)]))

