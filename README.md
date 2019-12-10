# R pipeline to predict causal variants underlying eQTLs

Author: Stephanie Hoyt
Date: July 2019

### Overview of scripts

**rank_qtl.R:** This script takes in data about QTL peaks (more information about specific inputs below) and annotates them with information about the clustering of strains at the peak, how many SNPs match that clustering pattern, how far the LOD peak is from the gene, what the gene biotype is, and whether or not the LOD peak and the gene start site are in the same TAD. This script does NOT return a score; it returns a dataframe of QTLs in the format expected by the scoring script.

**score_qtl.R:** This script takes in the dataframe produced by rank_qtl and assigns each QTL a numeric score based on the recorded annotations. The weights of each aspect of the score are passed into the function (lod_w, biallelic_w, distance_w, biotype_w, tad_w) and can be changed, or left at their default values. General idea behind scoring is that LOD score, and being bialleic will increase the score, the number of matching SNPs will increase the score up to a certain threshold, at which point it will start to decrease the score, and a greater distance from the peak to the gene will decrease the score. I’m also now including the gene_biotype information (being protein coding is positive and everything else is neutral).

**rank_snps.R:** The script takes in the dataframe produced by rank_qtl, as well as additional SNP annotation files (more information about specific inputs below), and annotates each of the SNPs found to match the pattern at each of the QTLs. This annotation includes assigning a chromHMM state, an ORegAnno state, checking whether the SNP is in the same TAD as the gene start site, getting the consequence from ensembl, calculating the distance from the SNP to the gene start site, and calculating the association between the SNP and the gene expression levels. This script does NOT return a score; it returns a list of QTLs and their SNPs in the format expected by the scoring script.

**score_snps.R:** This script takes in the dataframe produced by rank_snps.R and assigns each SNP a score based on the recorded annotations. Each of the chromHMM states, and each of the consequences are given a different score based on how much they indicate the SNP is likely to influence gene expression. The SNP association LOD score, the ORegAnno state, and the SNP being in the same TAD as the gene will increase the score. The chromHMM states and the consequences will also increase the score, but some more than others. The distance from the SNP to the gene start site will decrease the score.

### Input files and formats

###### For rank_qtl.R
All files should be passed in as character strings on command line, with the path coming first, followed by TAD file, then RData files (explained further below) and finally the name for the output RData file.

The TAD file originally used is from [the Yue lab](http://promoter.bx.psu.edu/hi-c/publications.html). TAD files from other sources can be used as input as long as they follow the same format.

These objects can be loaded in as part of RData files (multiple or all in one) named as command line arguments.
- covar: matrix object with the sample IDs as rownames and column indicating sex w/ 0 or 1 and any other covars as additional columns.
- exprZ: matrix object with normalized gene expression values
- kinship_loco: list generated with `qtl2::calc_kinship`
- probs: generated with `qtl2::calc_genoprob`
- gmap: genetic map
- peaks.annotated: dataframe required to have the following named columns in this precise order: "ensembl_gene_id", "peak_chr", "peak_cM", "lod", "ci_lo", "ci_hi", "interp_bp_peak", "before", "after", "gene_chrom", "gene_start", "gene_end", "strand", "mgi_symbol". Locations should be in bp, not Mbp.
- effects_blup: dataframe of estimated QTL effects with columns A-H for each of the 8 founder strains. output of `qtl2::scan1blup`

The file `cc_variants.sqlite` should be present in the same folder as other files, but does not need to be read in on the command line, because it is used for multiple cell types. Available [here](https://figshare.com/articles/SQLite_database_of_variants_in_Collaborative_Cross_founder_mouse_strains/5280229/2)

###### For score_qtl.R
First command line argument should be path, second will be the RData file output from `rank_qtl.R`, and the last should be the output file. No other input is required for this script because all annotations have been completed and are in the RData file.

###### For rank_snps.R
First command line argument should be the path, second should be TAD file name, third should be chromHMM file name (more details below), followed by RData files with specific objects (outlined below), and the last is the name of the output file.

TAD file should be the same as used in `rank_qtl.R`, and described above.

The chromHMM states file originally used is from [Guifeng Wei](https://github.com/guifengwei/ChromHMM_mESC_mm10). it was generated using the [ChromHMM program](http://compbio.mit.edu/ChromHMM/), so output from that program could be used as well, if the program is being applied to a different cell type.

These objects can be loaded in as part of RData files (multiple or all in one) named as command line arguments. Some are the same objects as passed into `rank_qtl.R` so will not be described again.
- See above for: covar, exprZ, kinship_loco, probs
- map_dat2: a dataframe with the columns (in this order) "chr", "pos", "marker", "chrom", "pos_bp", "n"
- peaks_local_over7: output from `rank_qtl.R`

The file `ORegAnno_Combined_2016.01.19.tsv` should be present in the same folder as other files, but does not need to be read in on the command line, because it is used for multiple cell types. Available [here](http://www.oreganno.org/dump/). The file `cc_variants.sqlite` is also still needed, as described above.


###### For score_snps.R
First command line argument should be path, next will be the RData file output from `rank_snps.R` and the RData file containing the possible consequences, and the last should be the output file. The consequences input is required so that each possible result can be annotated with a weight.


### Additional Documentation
- ranking.Rmd
- tads.Rmd
- final paper

### Examples

###### to run rank_qtl.R

`system("Rscript C:/Users/c-hoyts/Documents/GitHub/jax-ssp19/QTL_mapping/modular/scripts/rank_qtl.R C:/Users/c-hoyts/Documents/GitHub/jax-ssp19/QTL_mapping/modular/ mESC.Bonev_2017-raw.domains DO_mESC_paired_eQTL_forMapping.RData ENSMUSGid_to_symbol_v91.RData peaks_for_Stephanie.RData ESC_eQTL_effects.RData dummy_output.RData")`

###### to run score_qtl.R
`system("Rscript C:/Users/c-hoyts/Documents/GitHub/jax-ssp19/QTL_mapping/modular/scripts/score_qtl.R C:/Users/c-hoyts/Documents/GitHub/jax-ssp19/QTL_mapping/modular/ ranking_save_tads.RData scored_qtls.RData")`

###### to run rank_snps.R

`system("Rscript C:/Users/c-hoyts/Documents/GitHub/jax-ssp19/QTL_mapping/modular/scripts/rank_snps.R C:/Users/c-hoyts/Documents/GitHub/jax-ssp19/QTL_mapping/modular/ mESC.Bonev_2017-raw.domains mESC_E14_12_segments.bed.gz DO_mESC_paired_eQTL_forMapping.RData ranking_save_tads.RData test_score_snps.RData")`

###### to run score_snps.R

`system("Rscript C:/Users/c-hoyts/Documents/GitHub/jax-ssp19/QTL_mapping/modular/scripts/score_snps.R C:/Users/c-hoyts/Documents/GitHub/jax-ssp19/QTL_mapping/modular/ test_score_snps.RData possible_cons.RData scored_snps.RData")`
