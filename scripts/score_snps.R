library(tidyverse)

## Load Data
# args = commandArgs(trailingOnly = TRUE)
# path <- args[1]
# for(i in 2:length(args)-1){
#   load(paste0(path, args[i]))
# }
#path <- "/home/c-hoyts/qtl_data/modular/"
path <- "C:/Users/c-hoyts/Documents/GitHub/jax-ssp19/QTL_mapping/modular/"
load(paste0(path, "possible_cons.RData"))
load(paste0(path, "ranked_snps_w_oreganno_664.RData"))

########################
chromHMM_scores <- cbind(c("E1", "E2","E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12"),
                         c(1, 0, .2, 1, .5, 1, 1, 1, .5, .5, .8, 0)) %>% as.data.frame

#need to check on the 0
consequence_scores <- cbind(consequences, 
                            c(.7, .4, .7, 1, .7, .7, .4, .7, 1, 1, .4, .2, .4, .7, .4,
                              .7, .2, 1, .7, .4, .4, .4, .7, .4, .2, .7, .4, .2, .7, 
                              0, .4, .7, 1, 1, 1, .7, .7)) %>% as.data.frame()
colnames(consequence_scores)[2] <- "weight"


# score based on:
  # the consequence [positive, but some more than others]
  # oreganno state [positive]
  # chromHMM state [positive, but some more than others]
  # distance from SNP to gene [negative]
  # association LOD score [positive]
  # same tad as gene [positive]

score_snps <- function(snp, lod_w = .5, dist_w = 25, oreganno_w = 50, chromHMM_w = 25, tad_w = 15, conseq_w = 15, 
                       chromHMM_weights = chromHMM_scores, cons_weights = consequence_scores){
  snp <- snp %>% unlist
  if(is.na(snp["Assoc_LOD"])){
    lod <- 0
  }else{
    lod <- (snp["Assoc_LOD"] %>% as.character %>% as.numeric) * lod_w ## introduce a threshold? only count LOD > 7?
  }
  dist <- (snp["dist"] %>% as.character %>% as.numeric) * dist_w
  oreganno <- 0
  tad <- 0
  chromHMM <- 0
  conseq <- 0
  if(!(is.na(snp["oreganno_states"]))){
    oreganno <- oreganno_w
  }
  if(!(is.na(snp["snp_tad"]) | is.na(snp["gene_tad"]))){
    if(snp["snp_tad"] == snp["gene_tad"]){
      tad <- tad_w
    }
  }
  if(!(is.na(snp["chromHMM_states"]))){
    chromHMM <- (chromHMM_weights %>% filter(V1 == (snp["chromHMM_states"] %>% as.character)) %>% select(V2) %>% as.numeric) * chromHMM_w
  }
  if(!(is.na(snp["consequence"]))){
    cons <- str_split(snp["consequence"], ",") %>% unlist
    max <- 0
    for(i in 1:length(cons)){
      new <- (cons_weights %>% filter(SO_term == (cons[i] %>% as.character)) %>% select("weight") %>% as.numeric)
      if(new > max){
        max <- new
      }
    }
    conseq <- max * conseq_w
  }
  score <- (lod - dist + oreganno + tad + chromHMM + conseq)
  return(score)
}

score_all_snps <- function(all_snps){
  score_list <- all_snps %>% t %>% as.data.frame %>% map(score_snps)
  score <- score_list %>% unlist
  return(cbind(all_snps, score))
}

scores <- ranks %>% map(score_all_snps)

hmax_scores <- vector()
all_scores <- vector()
for(i in 1:length(scores)){
  all_scores <- c(all_scores, scores[[i]] %>% select(score) %>% unlist)
  max_scores <- c(max_scores, scores[[i]] %>% select(score) %>% unlist %>% max)
}

all_scores %>% hist(breaks = 50, main = "Distribution of SNP scores for eQTLs w/ fewer than 150 variants", xlim = c(-200,200), xlab = "Scores")

all_states <- vector()
for(i in 1:length(scores)){
  all_states <- c(all_states, scores[[i]] %>% select(oreganno_states) %>% unlist %>% unique)
}

all_states %>% unique()
#save(scores, chromHMM_scores, consequence_scores, file = paste0(path, "scored_snps.RData"))
