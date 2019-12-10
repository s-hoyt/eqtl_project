library(qtl2)
library(mclust)
library(parallel)
library(biomaRt)
library(tidyverse)

## Load Data
#read in data from command line instead of hard coding
#first argument should be path, second is TADs, then the next are RData files (order shouldn't matter). Last is output file
args = commandArgs(trailingOnly = TRUE)
#path <- "C:/Users/c-hoyts/Documents/GitHub/jax-ssp19/QTL_mapping/modular/"
path <- args[1]
tads <- read_tsv(paste0(path, args[2]), col_names = FALSE) %>% as.data.frame
for(i in 3:(length(args)-1)){
  load(paste0(path, args[i]))
}

## make functions to query SNPs and genes
#use hard coded string bc this should be constant across data sets
query_variants <- create_variant_query_func(paste0(args[1], "cc_variants.sqlite"))


#includes: covar, (expr), exprZ, kinship_loco, probs, (raw.expr), gmap, (map_dat2), (pmap), (covarTidy)
#load(paste0(path, "DO_mESC_paired_eQTL_forMapping.RData"))
#includes peaks.annotated
load(paste0(path,"peaks_for_Stephanie.RData"))
#includes effects_blup, (effects_std), (peaks), (effects.cmd)
#load(paste0(path, "ESC_eQTL_effects.RData"))


## Filtering

## add local/distant annotation
peaks_local_dis <- peaks.annotated %>% mutate(
  mid_gene = (gene_start+gene_end)/2, 
  loc_dis = ifelse((peak_chr != gene_chrom) | (abs(interp_bp_peak-mid_gene) > 10000000) , 'distant', 'local')) %>% cbind(effects_blup)

## add gene_biotype annotation
ensembl <- useMart(biomart="ensembl", host="http://dec2017.archive.ensembl.org","mmusculus_gene_ensembl")
gene_list <- peaks_local_dis %>% select(ensembl_gene_id) %>% distinct()
gene_list.w.bio <- getBM(attributes=c("ensembl_gene_id","external_gene_name","gene_biotype"),
                        filter='ensembl_gene_id', values=gene_list$ensembl_gene_id, mart=ensembl)

peaks_local_dis <- left_join(peaks_local_dis, gene_list.w.bio[,c(1,3)])

## just look at local peaks
## filter out only those LOD peaks over 7 and arrange from high to low
peaks_local_over7 <- peaks_local_dis %>% filter(loc_dis == "local", lod >= 7) %>% arrange(desc(lod))

## add column about distance from LOD peak to gene start using mutate
peaks_local_over7 <- peaks_local_over7 %>% mutate(peak2gene = abs(interp_bp_peak-gene_start))

## Clustering

#calculate clusters as given by mclust
mclust <- peaks_local_over7[, 17:24] %>% t() %>% as.data.frame() %>% 
  purrr::map(Mclust, G = c(2:4), verbose = FALSE)

clusters_mclust <- mclust %>% purrr::map(~ .$classification)

bic_mclust <- mclust %>% purrr::map(~ .$BIC)

## add clusters to peaks matrix using mutate
peaks_local_over7 <- peaks_local_over7 %>% mutate(mclust = clusters_mclust)

## add column indicating whether pattern is biallelic based on mclust results
peaks_local_over7 <- peaks_local_over7 %>% rowwise %>% mutate(biallelic = !(3 %in% unlist(mclust)))

#TADS

## get index column, name columns, and make chromosomes match other files
colnames(tads) <- c("chr", "start", "end")
tads <- tads %>% mutate(chr=replace(chr, chr == 23, "X"), chr = replace(chr, chr == 24, "Y"))


## Function to get matching TADs
find_tad <- function(snp, tads = tad_data){
  tads <- tads %>% filter(chr == snp[1])
  pos <- (snp[2] %>% as.numeric) * 1000000
  for(i in 1:nrow(tads)){
    if(between(pos, tads[i, 2], tads[i, 3])){
      return(i)
    }
  }
  return(NA)
}

## Call function and bind results
peak <- peaks_local_over7[,c("peak_chr", "interp_bp_peak")] %>% t %>% as.data.frame %>% map(find_tad, tads) %>% unlist
genes <- peaks_local_over7[,c("peak_chr","gene_start")] %>% t %>% as.data.frame %>% map(find_tad, tads) %>% unlist
print("past_tads")

peaks_local_over7 <- cbind(peaks_local_over7, peak, genes)


## Matching SNPs


## pseudo biallelic: check if peak fits pseudo biallelic criteria
## if so, return the modified pattern, and the index to ignore
pseudo_biallelic <- function(peaks){
  pattern <- peaks$mclust %>% as.data.frame() %>% t()
  # if there are at most 2 strains in the middle cluster
  if(sum(pattern == 2) <= 2){
    indices <- which(pattern == 2)
    avg <- peaks[indices+16] %>% unlist %>% mean %>% abs
    # see if they are close enough to zero to drop
    if(avg <= .1){
      # return modified cluster & which indices to ignore
      pattern[which(pattern == 2)] <- 0
      pattern[which(pattern == 3)] <- 2
      return(list(pattern, indices))
    }
  # if there are at most 2 strains in the 1st of 3rd cluster, merge them into the middle cluster
  }else if(sum(pattern == 1) <= 2 | sum(pattern == 3) <= 2){
    # return modified cluster and out of bounds index so as to still check all
    cluster <- case_when(sum(pattern == 1) <= 2 ~ 1,
                         sum(pattern == 3) <= 2 ~ 3)
    pattern[which(pattern == cluster)] <- 2
    pattern[which(pattern == 3)] <- 1
    return(list(pattern, 100))
  }
  # return NA for modified cluster to indicate no change, and that we shouldn't look for matching SNPs
  return(list(NA, NA))
}

## find_pattern: function to return the number of snps that match the pattern at a given peak
## returning lists instead of vectors so I have two discrete elements to work with
## checks for pseudo-biallelic by calling pseudo_biallelic()
## NOTE: no longer saving consequences here, because they won't line up nicely, but am saving the SNP position. 
## Can get consequences from ensembl, but can't get SNP position
find_pattern <- function(peaks){
  ## standard vars we need
  peak_loc <- peaks$interp_bp_peak/1000000
  ## default values of patterns & index:
  mod_pattern <- NA # indicates no change
  pattern <- peaks$mclust %>% as.data.frame() %>% t() # original clustering
  index <- 100 # out of range, so will not change vector later
  ## if the peak meets the pseudo biallelic parameters, set the index to skip and figure out which cluster is being dropped
  if(!(peaks$biallelic | 4 %in% peaks$mclust)){
    output <- pseudo_biallelic(peaks)
    if(!is.na(output[[2]])){
      pattern <- output[[1]] #used later
      mod_pattern <- output[[1]] # just here to be returned, signalling a change in pattern
      index <- output[[2]]
    }else{
      # pattern wasn't modified, so its still not biallelic and we won't be able to match SNPs
      return(list(-1, NA, NA, NA))
    }
  ## if it's just non-biallelic, then we won't be able to match any SNPs to it, so return
  }else if(!(peaks$biallelic)){
    return(list(-1, NA, NA, NA))
  }
  # make alt pattern (this is same for both biallelic and pseudo biallelic)
  # doing this because our cluster 1 is just the group with the lower effect, whereas cluster 1 in DB is reference allele
  # so we want to check both versions of the pattern, because we don't know which is reference
  ones <- which(pattern == 1)
  twos <- which(pattern == 2)
  alt_pattern <- vector(length = length(pattern))
  alt_pattern[ones] <- 2
  alt_pattern[twos] <- 1
  alleles <- query_variants(peaks$peak_chr, peak_loc-2, peak_loc+2)
  matches <- 0
  snps <- vector()
  position <- vector()
  ## loop through alleles to check for matches
  for(j in 1:nrow(alleles)){
    ## vector of indices to check
    check <- 8:15
    # if index is out of bounds (true biallelic) then the vector will not change and all strains will be checked
    check <- check[-index]
    if((alleles[j,check] == pattern) %>% all() | (alleles[j,check] == alt_pattern) %>% all()){
      matches <- matches + 1
      snps <- c(snps, alleles[j,1])
      position <- c(position, alleles[j,3])
    }
  }
  return(list(matches, position, snps, mod_pattern))
}

## pull out mod_pattern, or snps column
get_col <- function(result, num){
  column <- list()
  for(i in 1:length(result)){
    val <- result[[i]][num]
    if(is.null(val)){
      column[[i]] <- NA
    }else{
      column[[i]] <- result[[i]][num] %>% unlist()
    }
  }
  return(column)
}

#make a data frame of what was once a list
make_df <- function(result){
  num_snps <- vector(length = length(result))
  mod_pattern <- get_col(result, 4)
  snp_ids <- get_col(result, 3)
  position <- get_col(result, 2)
  for(i in 1:length(result)){
    num_snps[i] <- result[[i]][1] %>% as.integer()
  }
  value <- cbind(num_snps, position, snp_ids, mod_pattern)
  return(value)
}

## function to plot any peak given chromosome number and gene id, and optional annotation (ie scores, clusters)
plot_CC <- function(chrom, gene, annotation = NULL){
  out <- scan1(probs, exprZ[, gene, drop=FALSE],
                   kinship_loco, addcovar=covar)
  eff <- scan1blup(probs[, chrom], exprZ[, gene, drop=FALSE],
                       kinship_loco[[chrom]], addcovar=covar)
  peak <- find_peaks(out, gmap, threshold=7) %>% filter(chr == chrom) %>% .$pos
  par(mar = c(8, 4, 4, 2) + 0.1)
  plot_coefCC(eff, map = gmap[chrom], scan1_output=out, main = gene)
  ## add vertical line at peak
  abline(v = peak, col = "red", lty=2, lwd = 2)
  ## add info about scoring & clustering
  mtext(text = annotation, side = 1, outer = TRUE)
}

#plot_CC(2, "ENSMUSG00000027014")

## un-parallelized version
#output <- peaks_local_over7[1,] %>% t() %>% as.data.frame() %>% purrr::map(find_pattern) %>% make_df

## parallelized
#peaks_local_over7 %>% .[1:10,] %>% t() %>% as.data.frame() %>% mclapply(find_pattern, mc.cores = 4)


# apply find_pattern function then make results into a dataframe and save
snps_and_pos <- peaks_local_over7 %>% t() %>% as.data.frame() %>% mclapply(find_pattern, mc.cores = 16) %>% make_df

## Formatting output

peaks_local_over7 <- peaks_local_over7 %>% cbind(snps_and_pos)

qtl_score <- peaks_local_over7 %>%
  select(ensembl_gene_id, peak_chr, lod, mgi_symbol, peak2gene, biallelic, num_snps, mod_pattern, gene_biotype, peak, genes)

save(bic_mclust, peaks_local_over7, snps_and_pos, qtl_score, file = paste0(path, args[length(args)]))
