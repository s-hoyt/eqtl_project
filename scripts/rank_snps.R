library(qtl2)
library(parallel)
library(biomaRt)
library(tidyverse)

## Load Data
#first argument should be path, second is TADs, third is chromHMM, then the next are RData files (order shouldn't matter). Last is output file
args = commandArgs(trailingOnly = TRUE)
path <- args[1]
tad_data <- read_tsv(paste0(path, args[2]), col_names = FALSE) %>% as.data.frame
#tad_data <- read_tsv(paste0(path, "mESC.Bonev_2017-raw.domains"), col_names = FALSE) %>% as.data.frame
bed <- read.table(paste0(path, args[3]),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")q
#bed <- read.table(paste0(path, "mESC_E14_12_segments.bed.gz"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
for(i in 4:(length(args)-1)){
  load(paste0(path, args[i]))
}
#path <- "/home/c-hoyts/qtl_data/modular/"
path <- "C:/Users/c-hoyts/Documents/GitHub/jax-ssp19/QTL_mapping/modular/"
#includes: covar, (expr), exprZ, kinship_loco, probs, (raw.expr), (gmap), map_dat2, (pmap), (covarTidy)
load(paste0(path, "DO_mESC_paired_eQTL_forMapping.RData"))
#includes: (bic_mclust), peaks_local_over7, (qtl_score), (snps_and_pos)
load(paste0(path, "ranking_save_tads.RData"))
load(paste0(path, "scored_qtls.RData"))

## Get QTLs
## look at qtls with positive number of SNPs
qtls_w_snps <- peaks_local_over7 %>% filter(num_snps > 0, num_snps <= 150) %>% 
  select(ensembl_gene_id, peak_chr, interp_bp_peak, gene_start, mgi_symbol, gene_biotype, num_snps, snp_ids, position, peak, genes)


## Set up various databases/files I will be searching through later

#read in oreganno data and get just mouse data (most recent build)
oreganno <- read_tsv(paste0(path, "ORegAnno_Combined_2016.01.19.tsv")) %>% 
 filter(Species == "Mus musculus", Build == "mm10") %>% as.data.frame
colnames(oreganno)[8] <- "state" #why am I changing this if I'm also setting the other one. why not make the other match this?

#ChromHMM
## get start and end in terms of Mbp
bed <- bed %>% mutate(Start = V2/1000000, End = V3/1000000) %>% select(V1, Start, End, V4)
colnames(bed)[1] <- "Chr"
colnames(bed)[4] <- "state"

#TADs
## name columns, and make chromosomes match other files
tad_data <- cbind(tad_data, seq(from = 1, to = nrow(tad_data), by = 1))
colnames(tad_data) <- c("chr", "start", "end", "index")
tad_data <- tad_data %>% mutate(chr=replace(chr, chr == 23, "X"), chr = replace(chr, chr == 24, "Y"))


## loop through and get the state for each SNP
#commented out code was my attempt to not use a for loop, but I couldn't get it to work, so going to come back to that later
get_state <- function(snp, states){
  chr_name <- paste0("chr", snp[1])
  states <- states %>% filter(Chr == chr_name)
  for(i in 1:nrow(states)){
    if(between(snp[2] %>% as.character %>% as.numeric, states[i, "Start"], states[i, "End"])){
      return(states[i, "state"])
    }
  }
  return(NA)
}

## Function to get matching TADs
#a few differences from the above function so using both for now but will try to merge later
find_tad <- function(snp, tads = tad_data){
  tads <- tads %>% filter(chr == snp[1])
  pos <- (snp[2] %>% as.character %>% as.numeric) * 1000000
  for(i in 1:nrow(tads)){
    if(between(pos, tads[i, 2], tads[i, 3])){
      return(i)
    }
  }
  return(NA)
}

## make functions to query SNPs so that I can fill in missing consequences
query_variants <- create_variant_query_func(paste0(path, "cc_variants.sqlite"))


## change pmap and probs to make SNP association work
## code from Selcan
attr(probs, "is_x_chr") <- NULL

map_dat2 <- map_dat2 %>% mutate(Mbp = pos_bp/1000000)

split_map <- function (map, chr_names = NULL)
{
  map <- reorder_map_table(map, chr_names = chr_names)
  pos <- as.numeric(map[, 2])
  chr <- map[, 1]
  uchr <- unique(chr)
  names(pos) <- rownames(map)
  lapply(split(pos, factor(chr, uchr)), sort)
}

#This calls reorder_map_table() and the code for that is:
reorder_map_table <- function (map_tab, chr_col = 1, pos_col = 2, chr_names = NULL)
{
  chr <- map_tab[, chr_col]
  if (is.null(chr_names))
    chr_names <- unique(chr)
  chr <- factor(chr, levels = chr_names)
  pos <- suppressWarnings(as.numeric(map_tab[, pos_col]))
  map_tab[order(chr, pos, seq_len(nrow(map_tab))), , drop = FALSE]
}

pmap <- split_map(dplyr::select(map_dat2, marker,
                                chr, Mbp) %>% as.data.frame() %>%
                    tibble::remove_rownames() %>%
                    tibble::column_to_rownames('marker'))


## This is the main function
## Each QTL goes through this, and the SNPs for that QTL are all annotated, and then will be given a score 
#(function for score not yet done)
rank_snps <- function(qtl){
  ## get the important snp info
  snps_rank <- cbind(qtl$snp_ids %>% unlist, qtl$peak_chr %>% as.character, qtl$position %>% unlist) %>% as.data.frame()
  colnames(snps_rank) <- c("refsnp_id", "chr", "snp_pos")

  ## get consequences
  snps_cons <- vector(length = nrow(snps_rank))
  for(i in 1:nrow(snps_rank)){
    snps_cons[i] <- query_variants(snps_rank[i,2] %>% as.character, 
                                   snps_rank[i,3] %>% as.character %>% as.numeric - .001, 
                                   snps_rank[i,3] %>% as.character %>% as.numeric + .001) %>% 
      filter(snp_id == snps_rank[i,1]) %>% select(consequence)
  }
  snps_rank <- cbind(snps_rank, snps_cons %>% unlist)
  colnames(snps_rank)[4] <- "consequence"

  ## get rid of factors 
  #(not sure why they are factors in the first place, is there a way to stop this earlier?)
  snps_rank$chr <- snps_rank$chr %>% as.character
  snps_rank$snp_pos <- snps_rank$snp_pos %>% as.character %>% as.numeric
  
  ## add TAD classification to score (ie Same TAD, Adjacent, Distant to gene)
  snp_tad <- snps_rank[,c("chr", "snp_pos")] %>% t %>% as.data.frame %>% map(find_tad) %>% unlist
  gene_tad <- qtl$genes
  snps_rank <- cbind(snps_rank, snp_tad, gene_tad)

  ## add chromHMM states to the score
  chromHMM_states <- snps_rank[,2:3] %>% t %>% as.data.frame %>% map(get_state, states = bed) %>% unlist
  snps_rank <- cbind(snps_rank, chromHMM_states)

  ## ORegAnno
  #oreganno_states <- snps_rank[,2:3] %>% t %>% as.data.frame %>% map(get_state, states = oreganno) %>% unlist
  oreganno_states <- NA
  snps_rank <- cbind(snps_rank, oreganno_states)

  ## add distance from SNP to gene
  gene_pos <- qtl$gene_start/1000000
  snps_rank <- snps_rank %>% mutate(dist = abs(snp_pos - gene_pos))
  
  ## SNP association
  peak <- qtl$interp_bp_peak/1000000
  chr <- qtl$peak_chr
  
  snp_assoc <- scan1snps(probs, pmap, exprZ[,qtl$ensembl_gene_id], kinship_loco[[chr]], covar, 
                         query_func = query_variants, chr = chr, start = peak - 2, end = peak + 2, keep_all_snps = TRUE)$lod %>% 
    as.data.frame() %>% rownames_to_column()

  colnames(snp_assoc) <- c("refsnp_id", "Assoc_LOD")
  snps_rank <- left_join(snps_rank, snp_assoc)
  
  return(snps_rank)
}


#want a list of data frames
ranks <- qtls_w_snps %>% filter(ensembl_gene_id == "ENSMUSG00000048865") %>% t %>% as.data.frame %>% map(rank_snps)

#names(ranks) <- qtls_w_snps[c(3,9),1]

#save(ranks, file = paste0(path, args[length(args)]))