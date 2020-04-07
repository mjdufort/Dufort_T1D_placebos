##### scripts for analysis of data for Dufort et al. 2019. Cell type–specific immune phenotypes predict loss of insulin secretion in new-onset type 1 diabetes. (DOI: 10.1172/jci.insight.125556)
### this file includes scripts for differential expression analyses of individual genes using the limma package

##### set up environment: load packages #####

# load general packages
library(tidyverse)
library(magrittr)
theme_set(
  theme_bw(20) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour="black", fill=NA, size=1),
          axis.text=element_text(colour="black"),
          axis.ticks=element_line(colour="black")))
update_geom_defaults("point", list(shape=16))
library(gplots)
library(ordinal)

# load analysis_specific packages
library(edgeR)
library(limma)

# load custom packages (available at github.com/benaroyaresearch)
if (!require(countSubsetNorm)) remotes::install_github("benaroyaresearch/countSubsetNorm"); library(countSubsetNorm)
if (!require(RNAseQC)) remotes::install_github("benaroyaresearch/RNAseQC"); library(RNAseQC)
if (!require(limmaTools)) remotes::install_github("benaroyaresearch/limmaTools"); library(limmaTools)
if (!require(miscHelpers)) remotes::install_github("benaroyaresearch/miscHelpers"); library(miscHelpers)
if (!require(geneSetTools)) remotes::install_github("benaroyaresearch/geneSetTools"); library(geneSetTools)


##### load/save data from previous scripts #####

load("T1D_placebos_data_1_for_downstream_analyses.RData")
load("T1D_placebos_data_2_for_downstream_analyses.RData")


##### Set up and run limma for age_years #####

condition.tmp <- with(master.final, !is.na(age_years))
master.age_years <- fix_factors(filter(master.final, condition.tmp))
DGECounts.age_years.tmp <-
  calc_norm_counts(
    counts=counts.final, design=master.age_years, libID_col="libid",
    min_cpm=1, min_libs_perc=0.15,
    normalize=TRUE,
    return_DGEcounts=TRUE,
    group=master.age_years$participant_id)
master.age_years <-
  master.age_years[
    match(colnames(DGECounts.age_years.tmp),
          master.age_years[,"libid"]),]

DesignMat.age_years <-
  model.matrix(
    ~ age_years + sex + study_rna_batch,
    data=master.age_years)
vwts.age_years <-
  voomWithQualityWeights(
    DGECounts.age_years.tmp,
    design=DesignMat.age_years, plot=TRUE, span=0.1)
corfit.age_years.rand_patient <-
  duplicateCorrelation(
    vwts.age_years,
    design=DesignMat.age_years,
    block=master.age_years$participant_id)
vfit.age_years <-
  lmFit(vwts.age_years,
        block=master.age_years$participant_id,
        correlation=corfit.age_years.rand_patient$consensus.correlation) %>%
  eBayes()
topGenes.age_years <-
  topTable(vfit.age_years, coef = 2, number=Inf, sort.by="P")

## write out the full limma results
write.csv(
  topGenes.age_years %>%
    rownames_to_column(var="gene") %>%
    select(-threshold, -B),
  "Table_S3.csv",
  row.names=FALSE, quote=FALSE)

## volcano plot with B cell and neutrophil genes colored, based on Linsley gene modules
topGenes.tmp <- topGenes.age_years
gene_sets.tmp <- c("CD19.mod", "CXCR1.mod", "MPO.mod")
gene_set_colors.tmp <-
  c("red", "dodgerblue", "navyblue")
topGenes.tmp$linsley_gene_set <- NA
for (i in gene_sets.tmp) {
  topGenes.tmp$linsley_gene_set[
    rownames(topGenes.tmp) %in% gene_sets.Linsley[[i]]] <- i
}

## Volcano plot
plot_volcano_byvar_2var(
  topGenes.tmp[order(!is.na(topGenes.tmp$linsley_gene_set)),],
  point_order="input",
  plotdims=c(10,8),
  file_prefix="Fig_S8",
  x_lim=get_xy_lims(topGenes.tmp, x_symmetric=FALSE)[[1]],
  color_by_var="linsley_gene_set", color_var_lab="Gene set",
  color_by_var_levels=gene_sets.tmp,
  my_cols=gene_set_colors.tmp,
  p_cut=0.1, fc_cut=0)


##### Set up and run limma for cpeptide_slope #####

condition.tmp <-
  with(
    master.final,
    !is.na(log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept))
master.cpeptide_slope <-
  master.final %>%
  filter(condition.tmp) %>%
  fix_factors()
DGECounts.cpeptide_slope.tmp <-
  calc_norm_counts(
    counts=counts.final,
    design=master.cpeptide_slope,
    libID_col="libid",
    min_cpm=1, min_libs_perc=0.15, normalize=TRUE, return_DGEcounts=TRUE,
    group=master.cpeptide_slope$participant_id)
master.cpeptide_slope <-
  master.cpeptide_slope[
    match(colnames(DGECounts.cpeptide_slope.tmp),
          master.cpeptide_slope[,"libid"]),]

DesignMat.cpeptide_slope <-
  model.matrix(
    ~ log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept +
      sex + study_rna_batch,
    data=master.cpeptide_slope)
vwts.cpeptide_slope <-
  voomWithQualityWeights(
    DGECounts.cpeptide_slope.tmp,
    design=DesignMat.cpeptide_slope, plot=TRUE, span=0.1)
corfit.cpeptide_slope.rand_patient <-
  duplicateCorrelation(
    vwts.cpeptide_slope,
    design=DesignMat.cpeptide_slope,
    block=master.cpeptide_slope$participant_id)
vfit.cpeptide_slope <-
  lmFit(vwts.cpeptide_slope,
        block=master.cpeptide_slope$participant_id,
        correlation=corfit.cpeptide_slope.rand_patient$consensus.correlation) %>%
  eBayes()
topGenes.cpeptide_slope <-
  topTable(vfit.cpeptide_slope, coef = 2,
           number=Inf, sort.by="P")

# write out the full limma results
write.csv(
  topGenes.cpeptide_slope %>%
    rownames_to_column(var="gene") %>%
    select(-threshold, -B),
  "Table_S1.csv",
  row.names=FALSE, quote=FALSE)

## heatmap of genes related to DE processes
gene_sets.tmp <- c("CD19.mod", "MPO.mod", "CXCR1.mod")
p_cut.tmp <- 0.05; fc_cut.tmp <- 0
counts.tmp <-
  get_counts_sig_genes(
    vwts.cpeptide_slope,
    topGenes=topGenes.cpeptide_slope,
    p_cut=p_cut.tmp, fc_cut=fc_cut.tmp)
counts.tmp <-
  counts.tmp[
    rownames(counts.tmp) %in%
      unlist(gene_sets.Linsley[gene_sets.tmp]),]
plot_gene_heatmap(
  counts.tmp,
  design=master.cpeptide_slope,
  libID_col="libid",
  order_by_var="log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept",
  color_by_var="log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept",
  filename="Fig_S4.pdf",
  add_legend=TRUE, leg_x=0.92, leg_y=0.5,
  key=TRUE, keysize=0.8, density.info="none",
  labCol=FALSE)

## volcano plot with B cell and neutrophil genes colored, based on Linsley gene modules
topGenes.tmp <- topGenes.cpeptide_slope
gene_sets.tmp <- c("CD19.mod", "CXCR1.mod", "MPO.mod")
gene_set_colors.tmp <-
  c("red", "dodgerblue", "navyblue")
topGenes.tmp$linsley_gene_set <- NA
for (i in gene_sets.tmp) {
  topGenes.tmp$linsley_gene_set[
    rownames(topGenes.tmp) %in% gene_sets.Linsley[[i]]] <- i
}
plot_volcano_byvar_2var(
  topGenes.tmp[order(!is.na(topGenes.tmp$linsley_gene_set)),],
  point_order="input",
  plotdims=c(11,9),
  file_prefix="Fig_2",
  x_lim=c(-0.87, 0.87),
  y_lim=c(0,7.6),
  color_by_var="linsley_gene_set", color_var_lab="Gene set",
  color_by_var_levels=gene_sets.tmp,
  my_cols=gene_set_colors.tmp,
  p_cut=0.05, fc_cut=NULL)


##### Run limma for cpeptide_slope.age_years #####

DesignMat.cpeptide_slope.age_years <-
  model.matrix(
    ~ log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept + sex +
      study_rna_batch + age_years,
    data=master.cpeptide_slope) %>%
  drop.coef()
vwts.cpeptide_slope.age_years <-
  voomWithQualityWeights(
    DGECounts.cpeptide_slope.tmp,
    design=DesignMat.cpeptide_slope.age_years, plot=TRUE, span=0.1)
corfit.cpeptide_slope.age_years.rand_patient <-
  duplicateCorrelation(
    vwts.cpeptide_slope.age_years,
    design=DesignMat.cpeptide_slope.age_years,
    block=master.cpeptide_slope$participant_id)
vfit.cpeptide_slope.age_years <-
  lmFit(vwts.cpeptide_slope.age_years,
        block=master.cpeptide_slope$participant_id,
        correlation=corfit.cpeptide_slope.age_years.rand_patient$consensus.correlation) %>%
  eBayes()
topGenes.cpeptide_slope.age_years <-
  topTable(vfit.cpeptide_slope.age_years,
           coef = 2, number=Inf, sort.by="P")

## write out the full limma results
write.csv(
  topGenes.cpeptide_slope.age_years %>%
    rownames_to_column(var="gene") %>%
    select(-threshold, -B),
  "Table_S4.csv",
  row.names=FALSE, quote=FALSE)

## volcano plot with B cell and neutrophil genes colored, based on Linsley gene modules
topGenes.tmp <- topGenes.cpeptide_slope.age_years
gene_sets.tmp <- c("CD19.mod", "CXCR1.mod", "MPO.mod")
gene_set_colors.tmp <-
  c("red", "dodgerblue", "navyblue")
topGenes.tmp$linsley_gene_set <- NA
for (i in gene_sets.tmp) {
  topGenes.tmp$linsley_gene_set[
    rownames(topGenes.tmp) %in% gene_sets.Linsley[[i]]] <- i
}
plot_volcano_byvar_2var(
  topGenes.tmp[order(!is.na(topGenes.tmp$linsley_gene_set)),],
  point_order="input",
  plotdims=c(11,9),
  file_prefix="Fig_3E",
  x_lim=c(-0.9, 0.9),
  y_lim=c(0,4),
  color_by_var="linsley_gene_set", color_var_lab="Gene set",
  color_by_var_levels=gene_sets.tmp,
  my_cols=gene_set_colors.tmp,
  p_cut=0.05, fc_cut=NULL)


##### Set up and run limma for cpeptide_slope.neutrophils_lymphocytes #####

condition.tmp <-
  with(
    master_cbc_merged.filtered_cbcs,
    !is.na(log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept) &
      !is.na(neutrophils) & !is.na(lymphocytes))
master.cpeptide_slope.neutrophils_lymphocytes <-
  master_cbc_merged.filtered_cbcs %>%
  filter(condition.tmp) %>%
  fix_factors()
DGECounts.cpeptide_slope.neutrophils_lymphocytes.tmp <-
  calc_norm_counts(
    counts=counts.filtered_cbcs,
    design=master.cpeptide_slope.neutrophils_lymphocytes,
    libID_col="libid",
    min_cpm=1, min_libs_perc=0.15,
    normalize=TRUE,
    return_DGEcounts=TRUE,
    group=master.cpeptide_slope.neutrophils_lymphocytes$participant_id)
master.cpeptide_slope.neutrophils_lymphocytes <-
  master.cpeptide_slope.neutrophils_lymphocytes[
    match(colnames(DGECounts.cpeptide_slope.neutrophils_lymphocytes.tmp),
          master.cpeptide_slope.neutrophils_lymphocytes[,"libid"]),]

DesignMat.cpeptide_slope.neutrophils_lymphocytes <-
  model.matrix(
    ~ log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept +
      neutrophils + lymphocytes + sex + study_rna_batch,
    data=master.cpeptide_slope.neutrophils_lymphocytes)
vwts.cpeptide_slope.neutrophils_lymphocytes <-
  voomWithQualityWeights(
    DGECounts.cpeptide_slope.neutrophils_lymphocytes.tmp,
    design=DesignMat.cpeptide_slope.neutrophils_lymphocytes, plot=TRUE, span=0.1)
corfit.cpeptide_slope.neutrophils_lymphocytes.rand_patient <-
  duplicateCorrelation(
    vwts.cpeptide_slope.neutrophils_lymphocytes,
    design=DesignMat.cpeptide_slope.neutrophils_lymphocytes,
    block=master.cpeptide_slope.neutrophils_lymphocytes$participant_id)
vfit.cpeptide_slope.neutrophils_lymphocytes <-
  lmFit(vwts.cpeptide_slope.neutrophils_lymphocytes,
        block=master.cpeptide_slope.neutrophils_lymphocytes$participant_id,
        correlation=corfit.cpeptide_slope.neutrophils_lymphocytes.rand_patient$consensus.correlation) %>%
  eBayes()
topGenes.cpeptide_slope.neutrophils_lymphocytes <-
  topTable(vfit.cpeptide_slope.neutrophils_lymphocytes, coef = 2,
           number=Inf, sort.by="P")

## write out the full limma results
write.csv(
  topGenes.cpeptide_slope.neutrophils_lymphocytes %>%
    rownames_to_column(var="gene") %>%
    select(-threshold, -B),
  "Table_S2.csv",
  row.names=FALSE, quote=FALSE)

## volcano plot with B cell and neutrophil genes colored, based on Linsley gene modules
topGenes.tmp <- topGenes.cpeptide_slope.neutrophils_lymphocytes
p_cut.tmp <- 0.05
gene_sets.tmp <- c("CD19.mod", "CXCR1.mod", "MPO.mod")
gene_set_colors.tmp <-
  c("red", "dodgerblue", "navyblue")
topGenes.tmp$linsley_gene_set <- NA
for (i in gene_sets.tmp) {
  topGenes.tmp$linsley_gene_set[
    rownames(topGenes.tmp) %in% gene_sets.Linsley[[i]]] <- i
}

plot_volcano_byvar_2var(
  topGenes.tmp[order(!is.na(topGenes.tmp$linsley_gene_set)),],
  point_order="input",
  plotdims=c(11,9),
  file_prefix="Fig_3B",
  x_lim=c(-1.1, 1.1),
  y_lim=c(0,7),
  color_by_var="linsley_gene_set", color_var_lab="Gene set",
  color_by_var_levels=gene_sets.tmp,
  my_cols=gene_set_colors.tmp,
  p_cut=p_cut.tmp, fc_cut=NULL)


##### Set up and run limma for cpeptide_slope.age_years.interaction #####

condition.tmp <-
  with(
    master.final,
    !is.na(log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept) &
      !is.na(age_years))
master.cpeptide_slope.age_years.interaction <-
  master.final %>%
  filter(condition.tmp) %>%
  fix_factors()
DGECounts.cpeptide_slope.age_years.interaction.tmp <-
  calc_norm_counts(
    counts=counts.final,
    design=master.cpeptide_slope.age_years.interaction,
    libID_col="libid",
    min_cpm=1, min_libs_perc=0.15, normalize=TRUE, return_DGEcounts=TRUE,
    group=master.cpeptide_slope.age_years.interaction$participant_id)
master.cpeptide_slope.age_years.interaction <-
  master.cpeptide_slope.age_years.interaction[
    match(colnames(DGECounts.cpeptide_slope.age_years.interaction.tmp),
          master.cpeptide_slope.age_years.interaction[,"libid"]),]

DesignMat.cpeptide_slope.age_years.interaction <-
  model.matrix(
    ~ log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept + age_years +
      log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept:age_years +
      sex + study_rna_batch,
    data=master.cpeptide_slope.age_years.interaction)
vwts.cpeptide_slope.age_years.interaction <-
  voomWithQualityWeights(
    DGECounts.cpeptide_slope.age_years.interaction.tmp,
    design=DesignMat.cpeptide_slope.age_years.interaction, plot=TRUE, span=0.1)
corfit.cpeptide_slope.age_years.interaction.rand_patient <-
  duplicateCorrelation(
    vwts.cpeptide_slope.age_years.interaction,
    design=DesignMat.cpeptide_slope.age_years.interaction,
    block=master.cpeptide_slope.age_years.interaction$participant_id)
vfit.cpeptide_slope.age_years.interaction <-
  lmFit(vwts.cpeptide_slope.age_years.interaction,
        block=master.cpeptide_slope.age_years.interaction$participant_id,
        correlation=corfit.cpeptide_slope.age_years.interaction.rand_patient$consensus.correlation) %>%
  eBayes()
topGenes.cpeptide_slope.age_years.interaction <-
  topTable(vfit.cpeptide_slope.age_years.interaction,
           coef = "log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept:age_years",
           number=Inf, sort.by="P")

# write out the full limma results
write.csv(
  topGenes.cpeptide_slope.age_years.interaction %>%
    rownames_to_column(var="gene"),
  "Table_S5.csv",
  row.names=FALSE, quote=FALSE)

# get goana functional enrichment
goana.cpeptide_slope.age_years.interaction.P0.05 <-
  goana(vfit.cpeptide_slope.age_years.interaction,
        geneid=
          annotables::grch38$entrez[
            match(rownames(vfit.cpeptide_slope.age_years.interaction),
                  annotables::grch38$symbol)],
        coef=2,
        FDR=0.05) %>%
  mutate(Adj.P.Up=p.adjust(P.Up, method="BH"),
         Adj.P.Down=p.adjust(P.Down, method="BH"))
goana.cpeptide_slope.age_years.interaction.P0.05 %>%
  arrange(P.Up) %>%
  write.table(
    "cpeptide_slope.age_years.interaction.genes_P0.05.up.goana_results.txt",
    row.names=FALSE, quote=FALSE,
    sep="\t")
goana.cpeptide_slope.age_years.interaction.P0.05 %>%
  arrange(P.Down) %>%
  write.table(
    "cpeptide_slope.age_years.interaction.genes_P0.05.down.goana_results.txt",
    row.names=FALSE, quote=FALSE,
    sep="\t")


##### save objects for downstream use #####

save(file="T1D_placebos_data_3_for_downstream_analyses.RData",
     list=ls_grep("master|DesignMat|vfit|topGenes|corfit"))
