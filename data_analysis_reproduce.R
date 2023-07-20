
source("housman_test_cfDNA.R")


#Read in

library(methrix)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ComplexHeatmap)
library(RColorBrewer)

##download data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185307

meth_files <- list.files("public_data/", pattern = "hg38.bed.gz", full.names = T)
hg38 <- methrix::extract_CPGs(ref_genome = "BSgenome.Hsapiens.UCSC.hg38")

coldata <- read.delim("annotation_cfDNA_ref.txt", header = F)
colnames(coldata) <- c("Sample", "Category", "ID")
rownames(coldata) <- coldata$Sample

# to be able to order the plots by Category
coldata_order <- order(coldata$Category)


meth_cf <- methrix::read_bedgraphs(
files = meth_files,
ref_cpgs = hg38,
stranded = TRUE,
collapse_strands = TRUE, zero_based = F,
chr_idx = 1, start_idx = 2,
end_idx = 3, cov_idx = 10, beta_idx = 11, coldata = coldata)


#methylation atlas from:
#Moss, J., Magenheim, J., Neiman, D. et al.
#Comprehensive human cell-type methylation atlas reveals origins of circulating cell-free DNA in health and disease.
#Nat Commun 9, 5068 (2018). https://doi.org/10.1038/s41467-018-07466-6


#https://github.com/nloyfer/meth_atlas/blob/master/full_atlas.csv.gz
atlas <- read.csv("N:/data/reference_methylome_cfDNA/atlas_450k/atlas.csv", row.names=1)
atlas <- atlas[-1,]

### 450k manifest
### downloaded from: https://zwdzwd.github.io/InfiniumAnnotation
HM450.hg38.manifest <- read.delim("HM450.hg38.manifest.tsv.gz")
rownames(HM450.hg38.manifest) <- HM450.hg38.manifest$Probe_ID
HM450.hg38.manifest$CpG_beg <- HM450.hg38.manifest$CpG_beg + 1
colnames(HM450.hg38.manifest)[1:3]  <- c("chr", "start", "end")
HM450.hg38.manifest <- HM450.hg38.manifest[HM450.hg38.manifest$chr %in% paste0("chr", 1:22),]
common_cpg <- intersect(rownames(HM450.hg38.manifest), rownames(atlas))

HM450.hg38.manifest <- HM450.hg38.manifest[common_cpg,]
HM450 <- makeGRangesFromDataFrame(HM450.hg38.manifest)
HM450 <- HM450[!duplicated(HM450),]

###Create methrix from 450k reference atlas


common_cpg <- intersect(rownames(HM450.hg38.manifest), rownames(atlas))
atlas <- atlas[common_cpg,]

atlas_450k <- methrix:::create_methrix(beta_mat = atlas, cov_mat = ceiling(atlas), cpg_loci = HM450.hg38.manifest[common_cpg, c("chr", "start", "end")],
                           is_hdf5 = FALSE, genome_name = "hg38", col_data = data.frame(samples=colnames(atlas), row.names=colnames(atlas)), h5_dir = NULL,
                           ref_cpg_dt = hg38$cpgs, chrom_sizes = hg38$contig_lens, desc = NULL)
atlas_450k <- atlas_450k[names(HM450),]



## Selecting the most variable sites


results <- list()
overlapping <- list()
#downsampling
meth_cf <- remove_uncovered(meth_cf)
# to avoid selecting the very sparsely covered CpG sites
meth_cf <- coverage_filter(meth_cf, cov_thr = 1, min_samples = 2)




atlas_450k_most_variable <- methrix::order_by_sd(atlas_450k)[1:10000,]

res <- houseman_cell_type(m=meth_cf, reference_meth=atlas_450k_most_variable, included_regions=HM450[rownames(atlas_450k_most_variable)], pivot=NULL,
                              included_cell_types=NULL)

results[["il_atlas_most_variable_10000"]] <- round(res[[1]], 3)
overlapping[["il_atlas_most_variable_10000"]] <- res[[2]]



atlas_450k_most_variable <- methrix::order_by_sd(atlas_450k)[1:1000,]

res <- houseman_cell_type(m=meth_cf, reference_meth=atlas_450k_most_variable, included_regions=HM450[rownames(atlas_450k_most_variable)], pivot=NULL,
                         included_cell_types=NULL)

results[["il_atlas_most_variable_1000"]] <- round(res[[1]], 3)
overlapping[["il_atlas_most_variable_1000"]] <- res[[2]]

atlas_450k_most_variable <- methrix::order_by_sd(atlas_450k)[1:100,]



res <- houseman_cell_type(m=meth_cf, reference_meth=atlas_450k_most_variable, included_regions=HM450[rownames(atlas_450k_most_variable)], pivot=NULL,
                         included_cell_types=NULL)

results[["il_atlas_most_variable_100"]] <- round(res[[1]], 3)
overlapping[["il_atlas_most_variable_100"]] <- res[[2]]

#how does it look with randomly selected sites?
set.seed(563513)
sel <- sample(1:nrow(atlas_450k), 10000)
res <- houseman_cell_type(m=meth_cf, reference_meth=atlas_450k[sel,],  pivot=NULL, included_regions = HM450[sel,],
                         included_cell_types=NULL)

results[["il_atlas_random"]] <- round(res[[1]], 3)


#Use the pivot table to combine the epithelial cells

pivot_base <- read.csv("epithelial.csv")
pivot_base$destination[pivot_base$destination==""] <- "Other"
  pivot <- data.frame(matrix(NA, ncol=length(unique(pivot_base$destination)), nrow = nrow(pivot_base)))
  colnames(pivot) <- unique(pivot_base$destination)
  rownames(pivot) <- pivot_base$source

  for (i in seq_along(1:nrow(pivot))){
    pivot[i,]  <- ifelse(pivot_base$destination[i]==colnames(pivot), 1, 0)}


  rownames(pivot) <- gsub(" ", "_", rownames(pivot))
  rownames(pivot) <- gsub("-", ".", rownames(pivot))
  colnames(pivot) <- gsub(" ", "_", colnames(pivot))

  atlas_450k_most_variable <- methrix::order_by_sd(atlas_450k)[1:10000,]
  res <- houseman_cell_type(m=meth_cf, reference_meth=atlas_450k_most_variable, included_regions=HM450[rownames(atlas_450k_most_variable)], pivot=pivot,
                           included_cell_types=NULL)

  results[["il_atlas_most_variable_combined_cell_types_10000"]] <- round(res[[1]], 3)


  atlas_450k_most_variable <- methrix::order_by_sd(atlas_450k)[1:1000,]

  res <- houseman_cell_type(m=meth_cf, reference_meth=atlas_450k_most_variable, included_regions=HM450[rownames(atlas_450k_most_variable)], pivot=pivot,
                           included_cell_types=NULL)

  results[["il_atlas_most_variable_combined_cell_types_1000"]] <- round(res[[1]], 3)

  atlas_450k_most_variable <- methrix::order_by_sd(atlas_450k)[1:100,]


  res <- houseman_cell_type(m=meth_cf, reference_meth=atlas_450k_most_variable, included_regions=HM450[rownames(atlas_450k_most_variable)], pivot=pivot,
                           included_cell_types=NULL)

  results[["il_atlas_most_variable_combined_cell_types_100"]] <- round(res[[1]], 3)

##Plot heatmaps


  column_ha = HeatmapAnnotation(Type = coldata$Category[coldata_order],
                                col = list(Type=c(Control="gold","Lung adenocarcinoma"="darkgreen")))
library(circlize)
  col_fun = colorRamp2(c(0, 0.5, 1), c("darkblue", "white", "darkred"))


  h1 <- Heatmap(t(results[["il_atlas_most_variable_combined_cell_types_10000"]][coldata_order,]),
                cluster_rows = F, cluster_columns = F, top_annotation = column_ha, show_column_names = F,
                row_title  = "10000 sites", column_title = "Highly variable sites", col = col_fun)
  h2 <- Heatmap(t(results[["il_atlas_most_variable_combined_cell_types_1000"]][coldata_order,]), row_title = "1000 sites",
                cluster_rows = F, cluster_columns = F, show_column_names = F, col = col_fun)
  h3 <- Heatmap(t(results[["il_atlas_most_variable_combined_cell_types_100"]][coldata_order,]),  row_title = "100 sites",
                cluster_rows = F, cluster_columns = F, show_column_names = T, col = col_fun, column_labels = coldata$ID[coldata_order])
  pdf("N:/members/Reka/20_cfDNA_methylome/cell_type_cf/plots/most_variable.pdf")

  draw(h1  %v% h2  %v%h3)
dev.off()

cor(original_proportions$`Megalodon normal lung`, results[["il_atlas_most_variable_combined_cell_types_10000"]][,"Epithelial_cells"][coldata_order])
cor(original_proportions$`Megalodon normal lung`, results[["il_atlas_most_variable_combined_cell_types_1000"]][,"Epithelial_cells"][coldata_order])
cor(original_proportions$`Megalodon normal lung`, results[["il_atlas_most_variable_combined_cell_types_100"]][,"Epithelial_cells"][coldata_order])


  h1 <- Heatmap(t(results[["il_atlas_most_variable_10000"]][coldata_order,]),
                cluster_rows = F, cluster_columns = F, top_annotation = column_ha, show_column_names = F, column_title =  "10000 most variable sites")
  h2 <- Heatmap(t(results[["il_atlas_most_variable_1000"]][coldata_order,]), column_title = "1000 most variable sites",
                cluster_rows = F, cluster_columns = F, show_column_names = F)
  h3 <- Heatmap(t(results[["il_atlas_most_variable_100"]][coldata_order,]),  column_title = "10000 most variable sites",
                cluster_rows = F, cluster_columns = F, show_column_names = F)

  draw(h1  %v% h2  %v%h3)




#  atlas_450k_most_variable <- methrix::order_by_sd(atlas_450k)[1:100000,]
#  mv_gr <- get_matrix(atlas_450k_most_variable, in_granges = T, add_loci = T)
#  m_cf_gr <-  get_matrix(meth_cf_l[1:1000000,], in_granges = T, add_loci = T)
#subsetByOverlaps(mv_gr, m_cf_gr)
# using only highly covered
meth_cf <- methrix::mask_methrix(meth_cf, low_count = 2)
meth_cf <- remove_uncovered(meth_cf)
meth_cf <- coverage_filter(meth_cf, cov_thr = 1, min_samples = 3)

atlas_450k_most_variable <- methrix::order_by_sd(region_filter(atlas_450k, get_matrix(meth_cf, in_granges = T, add_loci = T)))[1:10000,]

res <- houseman_cell_type(m=meth_cf, reference_meth=atlas_450k_most_variable, included_regions=HM450[rownames(atlas_450k_most_variable)], pivot=pivot,
                         included_cell_types=NULL)

results[["il_atlas_most_variable_combined_cell_types_high_cov"]] <- round(res[[1]], 3)

