## scripts to load and analyze C-peptide data, and load RNAseq data for the T1D placebos project

##### set up environment: load packages #####

## load general packages
library(readxl)
library(tidyverse)
theme_set(
  theme_bw(20) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour="black", fill=NA, size=1),
          axis.text=element_text(colour="black"),
          axis.ticks=element_line(colour="black")))
update_geom_defaults("point", list(shape=16))
library(RColorBrewer)
library(ggthemes)

## load analysis-specific packages
library(limma)
library(edgeR)
library(annotables)
library(AER) # package for tobit models
library(caret) # package for model evaluation and cross-validation

# load custom packages (available at github.com/mjdufort)
library(RNAseQC)
library(countSubsetNorm)
library(miscHelpers)
library(geneSetTools)


##### import detailed study schedules, with visit numbers, names, corrected weeks, etc. #####

study_schedules <-
  read_xlsx(
    path="T1D_trials_info/T1D_trial_placebos_schedule_aggregated.xlsx",
    sheet=1) %>%
  standardize_dimnames()
for (
  i in colnames(study_schedules)[
    !(colnames(study_schedules) %in% c("study", "notes")) &
    !str_detect(colnames(study_schedules), "name|scheduled")]) {
  if (is.character(study_schedules[[i]]))
    study_schedules[[i]] <- study_schedules[[i]] %>%
      str_replace_all("<0", "NA") %>%
      as.numeric()
}


##### load patient and clinical data #####

## patient data (including age, sex, HLA)
patient_data <- read.csv("T1D_placebos_patient_data.csv")

## CBC data
cbc.merged <- read.csv("T1D_placebos_cbc_data_processed_merged.csv")

## C-peptide data
cpeptide_data.merged <- read.csv("T1D_placebos_cpeptide_data_processed_merged.csv")


##### load RNA-seq data #####

### read in counts data
counts <- read.csv("T1D_placebos_counts_merged.csv")

# Trim Lib ID's to lib#### (lib plus numbers)
colnames(counts) <-
  colnames(counts) %>%
  str_extract(pattern="lib[0-9]+") %>%
  make.unique(sep="_")

# Keep protein coding genes with HGNC symbols, and drop non-protein-coding genes
rownames(counts) <- counts[,1]
counts <- counts[, -1]
counts_to_aggregate.tmp <- counts
counts_to_aggregate.tmp$HGNC.symbol <-
  get_HGNC(rownames(counts_to_aggregate.tmp), type="protein_coding")

## sum counts for duplicated HGNC symbols, and drop rows with HGNC.symbols==NA
counts.aggregated <-
  aggregate(
    counts_to_aggregate.tmp[,1:(ncol(counts_to_aggregate.tmp)-1)],
    by=list(counts_to_aggregate.tmp$HGNC.symbol), sum)
rownames(counts.aggregated) <-
  counts.aggregated$Group.1
counts.aggregated <-
  counts.aggregated[
    ,colnames(counts.aggregated) != "Group.1"]
dim(counts.aggregated)


## read in library prep and annotation data
rnaseq_annotation.merged <- read.csv("T1D_placebos_rnaseq_annotation_merged.csv")

##### load RNA-seq library metrics #####

# read metrics files
metrics.merged <-
  read.csv("T1D_placebos_combined_metrics_merged.csv") %>%
  standardize_dimnames()

# standardize columns
for (i in 2:ncol(metrics.merged)) {
  if (!is.numeric(metrics.merged[[i]])) {
    if (sum(na.omit(str_detect(metrics.merged[[i]], "%"))) > 0) {
      metrics.merged[[i]] <- as.numeric(str_replace(metrics.merged[[i]], "%", "")) / 100
    } else metrics.merged[[i]] <- as.numeric(metrics.merged[[i]])
  }
}

# Trim Lib ID's to lib#### (lib plus 4 digit)
metrics.merged$libid <-
  metrics.merged$libid %>%
  str_extract("lib[0-9]+") %>%
  make.unique(sep="_")


##### Make quality cuts on RNAseq data #####

## merge sample annotation and metrics objects
rnaseq_annotation_metrics.merged <-
  merge(rnaseq_annotation.merged,
        metrics.merged,
        by="libid", all=FALSE)

# establish QC cuts by study
qc_cuts.by_study <-
  data.frame(
    study=c("ABATE", "START", "T1DAL", "TN02", "TN05", "TN09"),
    read_cut=c(1e6, 1e6, 1e6, 1.5e6, 5e6, 1.5e6),
    align_cut=c(0.9, 0.8, 0.8, 0.8, 0.8, 0.8),
    median_cv_cut=c(1, 1, 1, 0.7, 1, 0.8))

lib.bad.tmp <- character()
for (i in 1:nrow(qc_cuts.by_study))
  lib.bad.tmp <-
  c(lib.bad.tmp,
    rnaseq_annotation_metrics.merged$libid[
      rnaseq_annotation_metrics.merged$study==qc_cuts.by_study$study[i] &
        ((rnaseq_annotation_metrics.merged$pf_hq_aligned_reads < qc_cuts.by_study$read_cut[i] ) |
           (rnaseq_annotation_metrics.merged$mapped_reads_w_dups < qc_cuts.by_study$align_cut[i] ) |
           (rnaseq_annotation_metrics.merged$median_cv_coverage > qc_cuts.by_study$median_cv_cut[i] ))])
# 16 problematic libraries

## Make quality control cuts based on metrics
metrics.merged.qc <-
  metrics.merged[
    metrics.merged$libid %nin% lib.bad.tmp &
      metrics.merged$libid %in% rnaseq_annotation_metrics.merged$libid,]
nrow(metrics.merged.qc)  # 477 libraries (QC removed 16, or 3.2% of libraries)

# Remove from counts data libraries that fail QC cuts
counts.merged.aggregated.qc <-
  counts.aggregated[
    , colnames(counts.aggregated) %in% metrics.merged.qc$libid]
# 477 libraries

## Remove from sample annotation data libraries that fail QC cuts
rnaseq_annotation.merged.qc <-
  rnaseq_annotation.merged[
    rnaseq_annotation.merged$libid %in%
      colnames(counts.merged.aggregated.qc),]
# 477 libraries

## Force order of libraries
rnaseq_annotation.merged.qc <-
  rnaseq_annotation.merged.qc %>%
  arrange(libid)
counts.merged.aggregated.qc <-
  counts.merged.aggregated.qc[
    , match(rnaseq_annotation.merged.qc$libid, colnames(counts.merged.aggregated.qc))]


##### load combined object with patient data, RNAseq, C-peptide, CBC data merged #####

master <- read.csv("T1D_placebos_master_annotation.csv")


##### plot log AUC over time by subject #####

xlim.cpeptide_AUC_figures <-
  c(min(as.numeric(with(master, diagnosis_study_day)), na.rm=TRUE) + 40,
    max(master$cpeptide_study_day, na.rm=TRUE) - 50)

## plot AUC by individual vs cpeptide_study_day, with detection limit and diagnosis window, for manuscript
pdf("Fig_1A.pdf", w=8, h=5.3)
ggplot(
  data=master[!is.na(master$log_auc2hr),],
  mapping=aes(x=cpeptide_study_day, y=log_auc2hr, group=participant_id)) +
  geom_polygon(
    inherit.aes=FALSE,
    data=data.frame(
      cpeptide_study_day=rep(range(with(master, diagnosis_date - day_0), na.rm=TRUE), each=2),
      log_auc2hr=
        (range(master$log_auc2hr, na.rm=TRUE) + c(-0.5,0.5))[c(1,2,2,1)]),
    mapping=aes(x=cpeptide_study_day, y=log_auc2hr),
    fill="gray80") +
  geom_line(size=0.3) +
  scale_x_continuous(breaks=seq(0,1000, by=250)) +
  scale_y_continuous(breaks=log(10^(-2:0)), labels=10^(-2:0)) +
  labs(x="Days in study", y="C-peptide 2hr AUC\n(nmol / L / min)") +
  geom_hline(yintercept=log(low_lim_detect), linetype="dashed") +
  coord_cartesian(
    xlim=xlim.cpeptide_AUC_figures,
    ylim=range(master$log_auc2hr, na.rm=TRUE))
dev.off()


##### split patients into age tertiles #####

## generate age splits from patients in master with RNAseq libraries, then apply those splits to other patients
## split patients into even-sized groups by age_years
age_split.tmp <- unique(master[!is.na(master$libid), c("participant_id", "age_years")])
age_split.tmp$age_split3 <- cut_number(age_split.tmp$age_years, n=3)
master$age_split3 <-
  age_split.tmp$age_split3[
    match(master$participant_id, age_split.tmp$participant_id)]

## apply the cuts to patients in master but lacking RNAseq data
cuts.tmp <-
  levels(master$age_split3) %>%
  str_replace_all("\\[|\\]|\\(|\\)", "") %>%
  str_extract("[0-9\\.]+(?=,)") %>%
  as.numeric() %>%
  `[`(-1)
master$age_split3[is.na(master$age_split3)] <-
  levels(master$age_split3)[
    findInterval(master$age_years[is.na(master$age_split3)], cuts.tmp) + 1]

## apply age_split to patient_data.merged
patient_data.merged$age_split3 <-
  levels(master$age_split3)[
    findInterval(patient_data.merged$age_years, cuts.tmp) + 1] %>%
  factor(levels=levels(master$age_split3))


##### determine visits with c-peptide previously at detection threshold, to exclude them from model fitting #####

master$prev_at_cpep_detect_thresh <- as.logical(NA)
master <-
  master %>%
  arrange(participant_id, cpeptide_study_day)

for (i in unique(master$participant_id)) { # iterate over patient
  rows.tmp <- which(master$participant_id==i & !is.na(master$auc2hr)) # determine all rows for current patient
  master$prev_at_cpep_detect_thresh[rows.tmp[1]] <- FALSE # set first visit to FALSE
  if (length(rows.tmp) > 1)
    for (j in (2:length(rows.tmp))) {
      master$prev_at_cpep_detect_thresh[rows.tmp[j]] <-
        (master$auc2hr[rows.tmp[j-1]] < low_lim_detect * 1.01) |  # previous visit is at detection limit (lowest value about detection limit is ~1.08 of detection limit)
        master$prev_at_cpep_detect_thresh[rows.tmp[j-1]]  # previous visit had prev_at_cpep_detect_thresh
    }
}

##### scale cpeptide_study_day so that it works better in linear models #####
cpeptide_study_day_scaled.tmp <-
  scale(master$cpeptide_study_day, center=FALSE,
        scale=sd(master$cpeptide_study_day, na.rm=TRUE))
master$cpeptide_study_day_scaled <- as.vector(cpeptide_study_day_scaled.tmp)
days.sd <- attr(cpeptide_study_day_scaled.tmp, "scaled:scale")


##### calculate rates of log C-peptide AUC change by cpeptide_study_day using linear models #####

## fit model to log_auc2hr that includes a random patient-specific intercept and slope

# fit models with first-order polynomial
lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random <-
  lmer(log_auc2hr ~ 1 + (1|participant_id) + cpeptide_study_day_scaled +
         (cpeptide_study_day_scaled-1|participant_id),
       data=
         master[
           with(master,
                !prev_at_cpep_detect_thresh & !duplicated(cpeptide_visit_id)),])


## plot model fits to C-peptide log AUC over time by subject
pdf("Fig_S1B.pdf", w=7.9, h=6)
par(mar=c(5,6.5,4,2)+0.1)
plot(x=0, type="n",
     xlim= (xlim.cpeptide_AUC_figures - c(10,0)) / days.sd,
     ylim=range(master$log_auc2hr, na.rm=TRUE),
     xlab="Days in study", ylab="C-peptide 2hr AUC\n(nmol / L / min)", xaxt="n", yaxt="n",
     cex.lab=1.5)
axis(1,
     at = seq(0,1000, by=250)/days.sd,
     labels = seq(0,1000, by=250),
     cex.axis=1.5)
axis(2, at = log(10^(-2:0)), labels=10^(-2:0),
     cex.axis=1.5, las=1)
polygon(
  x = rep(range(with(master, diagnosis_date - day_0), na.rm=TRUE), each=2) / days.sd,
  y = (range(master$log_auc2hr, na.rm=TRUE) + c(-0.5,0.5))[c(1,2,2,1)],
  col="gray80", border=NA)
for (i in unique(master$participant_id)) {
  if (
    sum(!is.na(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id[
      match(i, rownames(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id)),]))==
    ncol(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id)) {
    b0.tmp <-
      coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id$'(Intercept)'[
        match(i, rownames(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id))]
    b1.tmp <-
      coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id$cpeptide_study_day_scaled[
        match(i, rownames(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id))]
    plot(
      function(x) {b0.tmp + (x*b1.tmp)},
      from = min(master$cpeptide_study_day[master$participant_id==i], na.rm=TRUE)/days.sd, 
      to = max(master$cpeptide_study_day[master$participant_id==i], na.rm=TRUE)/days.sd,
      lwd=1,
      add=TRUE)
  }
}
abline(h=log(low_lim_detect), lty=2)
dev.off()


##### extract rates of C-peptide change from model #####

## extract slope at points based on the model
master$log_auc2hr_intercept_cpeptide_study_day_linear_random_with_intercept <- as.numeric(NA)
master$log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept <- as.numeric(NA)

for (i in 1:nrow(master)) {
  cat("Starting row", i, "of", nrow(master), "\n")
  
  x.tmp <- master$rnaseq_study_day[i] / days.sd # get the x-value for the data point
  
  # extract slope and intercept for linear random with intercept
  if (sum(!is.na(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id[
    match(
      master$participant_id[i],
      rownames(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id)),])) == 2) {
    b0.tmp <-
      coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id[["(Intercept)"]][
        match(master$participant_id[i], rownames(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id))]
    b1.tmp <-
      coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id[["cpeptide_study_day_scaled"]][
        match(master$participant_id[i], rownames(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id))]
    
    master$log_auc2hr_intercept_cpeptide_study_day_linear_random_with_intercept[i] <- b0.tmp
    master$log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept[i] <-
      b1.tmp / days.sd * 365.25
  }
}


##### plot C-peptide rate vs. age by study #####

## plot C-peptide rate vs. age for each study, with gray dots for other other studies, with smooth for each study
for (i in sort(unique(master$study))) {
  name.tmp <-
    switch(
      i,
      "ABATE"="Teplizumab",
      "START"="ATG",
      "T1DAL"="Alefacept",
      "TN02"="MMF / Daclizumab",
      "TN05"="Rituximab",
      "TN09"="Abatacept")
  pdf(
    paste0("Fig_S3_", i, ".pdf"),
    w=9, h=6)
  print(
    ggplot(
      mapping=aes(
        x=age_years, y=log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept)) +
      geom_point(
        data=master[!duplicated(master$participant_id) & master$study!=i,],
        color="gray60",
        size=2) +
      lims(
        x=range(master$age_years, na.rm=TRUE),
        y=range(master$log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept, na.rm=TRUE)) +
      geom_smooth(
        # data=master[!duplicated(master$participant_id),],
        data=master[!duplicated(master$participant_id) & master$study==i,],
        color="black", se=FALSE, linetype="dashed",
        method="lm", formula= y~I(log(x)),
        size=1) +
      geom_point(
        data=master[!duplicated(master$participant_id) & master$study==i,],
        color="black",
        size=5) +
      # labs(title=name.tmp, x="Age at diagnosis", y="Rate of C-peptide change"))
      labs(title=NULL, x=NULL, y=NULL) +
      geom_text(
        mapping=aes(label=study),
        data=data.frame(
          age_years=45,
          log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept=-3,
          study=name.tmp),
        size=12, hjust=1))
  dev.off()
}

##### split subjects into fast and slow progressors #####

cpep_rate.RNAseq_patients.tmp <-
  unique(
    master[
      !is.na(master$libid),
      c("participant_id",
        "log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept")])
any(duplicated(cpep_rate.RNAseq_patients.tmp$participant_id)) # none have multiple values of slope - good!

## fastest 25 vs. slowest 75 percent
master$progressor_split_cpep_rate_fast25_slow75 <- 
  ifelse(
    master$log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept >
      quantile(cpep_rate.RNAseq_patients.tmp$log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept,
               prob=0.25, na.rm=TRUE),
    "slow", "fast") %>%
  factor(levels=c("fast", "slow"))


# plot C-peptide log AUC by subject over time, with lines colored by fast vs. slow progressor
pdf("Fig_S1D.pdf", w=7.6, h=5.3)
ggplot(
  data=master[!is.na(master$log_auc2hr),],
  mapping=aes(
    x=cpeptide_study_day, y=log_auc2hr,
    group=participant_id,
    color=progressor_split_cpep_rate_fast25_slow75)) +
  geom_polygon(
    inherit.aes=FALSE,
    data=data.frame(
      cpeptide_study_day=rep(range(with(master, diagnosis_date - day_0), na.rm=TRUE), each=2),
      log_auc2hr=
        (range(master$log_auc2hr, na.rm=TRUE) + c(-0.5,0.5))[c(1,2,2,1)]),
    mapping=aes(x=cpeptide_study_day, y=log_auc2hr),
    fill="gray80") +
  geom_line(size=0.3) +
  scale_color_manual(
    "Progressor",
    values=c("fast"="red", "slow"="blue")) +
  guides(color=FALSE) +
  scale_x_continuous(breaks=seq(0,1000, by=250)) +
  scale_y_continuous(breaks=log(10^(-2:0)), labels=10^(-2:0)) +
  labs(x="Days in study", y="C-peptide 2hr AUC\n(nmol / L / min)") +
  geom_hline(yintercept=log(low_lim_detect), linetype="dashed") +
  coord_cartesian(
    xlim=xlim.cpeptide_AUC_figures,
    ylim=range(master$log_auc2hr, na.rm=TRUE))
dev.off()


##### plot histogram of C-peptide rate by progressor status #####

pdf("Fig_S1C.pdf", w=8, h=6)
ggplot(
  master %>%
    filter(!is.na(progressor_split_cpep_rate_fast25_slow75)) %>%
    select(participant_id, log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept,
           progressor_split_cpep_rate_fast25_slow75) %>%
    unique(),
  mapping=aes(
    x=log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept,
    group=progressor_split_cpep_rate_fast25_slow75, fill=progressor_split_cpep_rate_fast25_slow75)) +
  geom_histogram(position="identity", color="black") +
  scale_fill_manual("Progressor", values=c("red", "blue")) +
  labs(x="Rate of C-peptide change") +
  guides(fill=FALSE)
dev.off()


##### filter master and counts to include only shared libraries, and put them in the same order #####

# now filter the master and counts object, and put them in the same order
master.shared <-
  master[master$libid %in% colnames(counts.merged.aggregated.qc),] %>%
  arrange(participant_id, cpeptide_study_day, libid)

# drop duplicated samples
for (i in unique(master.shared$rnaseq_visit_id[duplicated(master.shared$rnaseq_visit_id)])) {
  libs.tmp <-
    master.shared$libid[
      master.shared$rnaseq_visit_id %in% i]
  lib_to_keep.tmp <-  # keep library with higher alignment rate
    libs.tmp[
      which.max(
        metrics.merged$mapped_reads_w_dups[
          match(libs.tmp, metrics.merged$libid)])]
  master.shared <-
    master.shared[
      (master.shared$rnaseq_visit_id %nin% i) |
        (master.shared$libid %in% lib_to_keep.tmp),]
}
duplicated(master.shared$rnaseq_visit_id)  # looks good

counts.shared <-
  counts.merged.aggregated.qc[
    ,match(master.shared$libid, colnames(counts.merged.aggregated.qc))]

# drop excess "scaled_" variables from master
master.shared <-
  master.shared[, !str_detect(colnames(master.shared), "scaled_")]

# force same library order
counts.shared <-
  counts.shared[,match(master.shared$libid, colnames(counts.shared))]


##### filter counts and run PCA on (normalized) full data set to check for bias and batch effects #####

counts.shared.filtered_genes <-
  countSubsetNorm::calc_norm_counts(
    counts=counts.shared, design=master.shared,
    libID_col="libid",
    min_cpm=1, min_libs_perc=0.15,
    normalize=FALSE,
    return_DGEcounts=FALSE)
# filtered to 12,877 genes

# run PCA on the filtered, non-normalized data
pcaAll.counts.shared.filtered_genes <-
  RNAseQC::calc_PCAs(
    as.data.frame(counts.shared.filtered_genes),
    log2_transform=TRUE)

# merge sample info with PC scores  
scores_design_pcaAll.counts.shared.filtered_genes <-
  merge(master.shared,
        as.data.frame(pcaAll.counts.shared.filtered_genes$x),
        by.x = "libid", by.y="row.names")


### Color plots of PCs by different variables

## plots of PC1s 1-3, colored by study
plot_PCAs(
  scores_design_pcaAll.counts.shared.filtered_genes,
  pvars.labs=pcaAll.counts.shared.filtered_genes$pvars.labs, PCs = 1:3,
  color_by_var="study",
  color_by_var_levels=sort(unique(scores_design_pcaAll.counts.shared.filtered_genes$study)),
  my_cols=colorblind_pal()(8)[c(1:4,6,8)],
  plotdims=c(9,6))

## plots of PC1s 1-3, colored by rnaseq batch variables, for each study
rna_batch_vars.tmp <- c("source", "rna_processing_core", "seq_site", "rna_batch")
for (i in unique(master.shared$study)) {
  scores_design_pcaAll.counts.shared.filtered_genes.tmp <-
    dplyr::filter(scores_design_pcaAll.counts.shared.filtered_genes, study==i)
  plot_vars.tmp <-
    rna_batch_vars.tmp[
      colSums(!is.na(scores_design_pcaAll.counts.shared.filtered_genes.tmp[, rna_batch_vars.tmp])) > 0]
  for (j in plot_vars.tmp)
    plot_PCAs(
      scores_design_pcaAll.counts.shared.filtered_genes.tmp,
      pvars.labs=pcaAll.counts.shared.filtered_genes$pvars.labs, PCs=1:3,
      color_by_var=j,
      my_cols=colorblind_pal()(8),
      plotdims=c(9,6))
}

## quantify correlations of annotation variables with PCs, by study, and output heatmaps
cor_all_clinvars_vs_pcaAll.counts.shared.filtered_genes.by_study <- list()
for (i in sort(unique(master.shared$study))) {
  cor_all_clinvars_vs_pcaAll.counts.shared.filtered_genes.by_study[[i]] <-
    calc_PCcors(pcaAll.counts.shared.filtered_genes, master.shared[master.shared$study==i,])
  
  try(plot_PCcor_heatmap(
    cor_all_clinvars_vs_pcaAll.counts.shared.filtered_genes.by_study[[i]][1:10,],
    plotdims=c(12,6)))
}
# in AbATE, PC2 and PC3 correlated with seq_site, a little with rna_processing_core (use rna_batch)
# in TN02, PC1 strongly correlated with all batch varaiables (use rna_batch)
# in TN09, PC1 is strongly correlated with both batch variables (use rna_batch)


##### generate batch variables for use in models #####

## library prep and sequencing batch stuff
# create empty variables for different prep step batches
# concatenate different prep step batches
# create study_batch variable to concatenate the study name (to start) plus all prep batch variables (to control for all this)
# merge all prep step batch data into master.shared for the studies where I have it
# replace study_batch variable with study where it's NA
# compare within each study by the different prep step batches

# for all studies except ABATE, TN02, and TN09, use study
# for ABATE, TN02, and TN09, use study_rna_batch

master.shared$study_rna_batch <-
  master.shared$study
master.shared$study_rna_batch[master.shared$study %in% c("ABATE", "TN02", "TN09")] <-
  with(master.shared[master.shared$study %in% c("ABATE", "TN02", "TN09"),],
       paste(study, rna_batch, sep="_"))


##### generate simply-named combined objects with all libraries passing QC and outlier cuts #####

## these objects have had bad libraries and outliers removed
## the data were normalized previously
counts.final <- counts.shared
master.final <- master.shared


##### generate vwts.all #####

## generate vwts.all object
vwts.all <-
  calc_norm_counts(
    counts=counts.final, design=master.final, libID_col="libid",
    min_cpm=1, min_libs_perc=0.15,
    normalize=TRUE,
    return_DGEcounts=TRUE) %>%
  voom(plot=TRUE, span=0.1)


##### remove batch effects from ABATE, TN02, and TN09, and check if that removes the batch effects #####

### for ABATE
# remove seq_site or rna_batch

vwts.ABATE <-
  vwts.all[, match(master.final$libid[master.final$study=="ABATE"],
                   colnames(vwts.all$E))]

vwts.ABATE.seq_site_removed <-
  removeBatchEffect(
    vwts.ABATE,
    batch=master.final$seq_site[master.final$study=="ABATE"])
pca.ABATE.seq_site_removed <-
  calc_PCAs(vwts.ABATE.seq_site_removed, cpm=FALSE, log2_transform=FALSE)
scores_design_pca.ABATE.seq_site_removed <-
  merge(
    master.final[master.final$study=="ABATE",],
    as.data.frame(
      pca.ABATE.seq_site_removed$x),
    by.x = "libid", by.y="row.names")

vwts.ABATE.rna_batch_removed <-
  removeBatchEffect(
    vwts.ABATE,
    batch=master.final$rna_batch[master.final$study=="ABATE"])
pca.ABATE.rna_batch_removed <-
  calc_PCAs(vwts.ABATE.rna_batch_removed, cpm=FALSE, log2_transform=FALSE)
scores_design_pca.ABATE.rna_batch_removed <-
  merge(
    master.final[master.final$study=="ABATE",],
    as.data.frame(
      pca.ABATE.rna_batch_removed$x),
    by.x = "libid", by.y="row.names")

# plot PCAs after batch effect removal, colored by batch variables
for (j in c("seq_site", "rna_processing_core", "rna_batch"))
  plot_PCAs(
    scores_design_pca.ABATE.seq_site_removed,
    pvars.labs=pca.ABATE.seq_site_removed$pvars.labs, PCs=1:3,
    color_by_var=j,
    my_cols=colorblind_pal()(8),
    file_prefix="plots/T1D_placebos_ABATE.seq_site_removed",
    plotdims=c(9,6))

for (j in c("seq_site", "rna_processing_core", "rna_batch"))
  plot_PCAs(
    scores_design_pca.ABATE.rna_batch_removed,
    pvars.labs=pca.ABATE.rna_batch_removed$pvars.labs, PCs=1:3,
    color_by_var=j,
    my_cols=colorblind_pal()(8),
    file_prefix="plots/T1D_placebos_ABATE.rna_batch_removed",
    plotdims=c(9,6))

# use pcCors to check removal of batch effect
cor_all_clinvars_vs_pcaAll.counts.ABATE.seq_site_removed <-
  calc_PCcors(
    pca.ABATE.seq_site_removed,
    master.final[master.final$study=="ABATE",])
try(plot_PCcor_heatmap(
  cor_all_clinvars_vs_pcaAll.counts.ABATE.seq_site_removed[1:10,],
  filename=
    "plots/T1D_plac.heatmap_correlations_clinvars_pcaAll.counts.ABATE.seq_site_removed.pdf",
  plotdims=c(12,6)))

cor_all_clinvars_vs_pcaAll.counts.ABATE.rna_batch_removed <-
  calc_PCcors(
    pca.ABATE.rna_batch_removed,
    master.final[master.final$study=="ABATE",])
try(plot_PCcor_heatmap(
  cor_all_clinvars_vs_pcaAll.counts.ABATE.rna_batch_removed[1:10,],
  filename=
    "plots/T1D_plac.heatmap_correlations_clinvars_pcaAll.counts.ABATE.rna_batch_removed.pdf",
  plotdims=c(12,6)))

## for ABATE, rna_batch does much better, especially with eliminating a residual batch effect in rna_processing_core


### for TN02
## remove source and rna_processing_core, or rna_batch

vwts.TN02 <-
  vwts.all[
    , match(master.final$libid[master.final$study=="TN02"],
            colnames(vwts.all$E))]

vwts.TN02.source_and_rna_processing_core_removed <-
  removeBatchEffect(
    vwts.TN02,
    batch=master.final$source[master.final$study=="TN02"],
    batch2=master.final$rna_processing_core[master.final$study=="TN02"])
pca.TN02.source_and_rna_processing_core_removed <-
  calc_PCAs(vwts.TN02.source_and_rna_processing_core_removed,
            cpm=FALSE, log2_transform=FALSE)
scores_design_pca.TN02.source_and_rna_processing_core_removed <-
  merge(
    master.final[master.final$study=="TN02",],
    as.data.frame(
      pca.TN02.source_and_rna_processing_core_removed$x),
    by.x = "libid", by.y="row.names")

vwts.TN02.rna_batch_removed <-
  removeBatchEffect(
    vwts.TN02,
    batch=master.final$rna_batch[master.final$study=="TN02"])
pca.TN02.rna_batch_removed <-
  calc_PCAs(vwts.TN02.rna_batch_removed, cpm=FALSE, log2_transform=FALSE)
scores_design_pca.TN02.rna_batch_removed <-
  merge(
    master.final[master.final$study=="TN02",],
    as.data.frame(
      pca.TN02.rna_batch_removed$x),
    by.x = "libid", by.y="row.names")

# now plot them
for (j in c("source", "rna_processing_core", "rna_batch"))
  plot_PCAs(
    scores_design_pca.TN02.source_and_rna_processing_core_removed,
    pvars.labs=pca.TN02.source_and_rna_processing_core_removed$pvars.labs, PCs=1:3,
    color_by_var=j,
    my_cols=colorblind_pal()(8),
    file_prefix="plots/T1D_placebos_TN02.source_and_rna_processing_core_removed",
    plotdims=c(9,6))

for (j in c("source", "rna_processing_core", "rna_batch"))
  plot_PCAs(
    scores_design_pca.TN02.rna_batch_removed,
    pvars.labs=pca.TN02.rna_batch_removed$pvars.labs, PCs=1:3,
    color_by_var=j,
    my_cols=colorblind_pal()(8),
    file_prefix="plots/T1D_placebos_TN02.rna_batch_removed",
    plotdims=c(9,6))

# use pcCors to check removal of batch effect
cor_all_clinvars_vs_pcaAll.counts.TN02.source_and_rna_processing_core_removed <-
  calc_PCcors(
    pca.TN02.source_and_rna_processing_core_removed,
    master.final[master.final$study=="TN02",])
try(plot_PCcor_heatmap(
  cor_all_clinvars_vs_pcaAll.counts.TN02.source_and_rna_processing_core_removed[1:10,],
  filename=
    "plots/T1D_plac.heatmap_correlations_clinvars_pcaAll.counts.TN02.source_and_rna_processing_core_removed.pdf",
  plotdims=c(12,6)))

cor_all_clinvars_vs_pcaAll.counts.TN02.rna_batch_removed <-
  calc_PCcors(
    pca.TN02.rna_batch_removed,
    master.final[master.final$study=="TN02",])
try(plot_PCcor_heatmap(
  cor_all_clinvars_vs_pcaAll.counts.TN02.rna_batch_removed[1:10,],
  filename=
    "plots/T1D_plac.heatmap_correlations_clinvars_pcaAll.counts.TN02.rna_batch_removed.pdf",
  plotdims=c(12,6)))

## for TN02, both work similarly well; use rna_batch


### for TN09
# remove source and rna_processing_core, or rna_batch

vwts.TN09 <-
  vwts.all[
    , match(master.final$libid[master.final$study=="TN09"],
            colnames(vwts.all$E))]

vwts.TN09.source_and_rna_processing_core_removed <-
  removeBatchEffect(
    vwts.TN09,
    batch=master.final$source[master.final$study=="TN09"],
    batch2=master.final$rna_processing_core[master.final$study=="TN09"])
pca.TN09.source_and_rna_processing_core_removed <-
  calc_PCAs(vwts.TN09.source_and_rna_processing_core_removed,
            cpm=FALSE, log2_transform=FALSE)
scores_design_pca.TN09.source_and_rna_processing_core_removed <-
  merge(
    master.final[master.final$study=="TN09",],
    as.data.frame(
      pca.TN09.source_and_rna_processing_core_removed$x),
    by.x = "libid", by.y="row.names")

vwts.TN09.rna_batch_removed <-
  removeBatchEffect(
    vwts.TN09,
    batch=master.final$rna_batch[master.final$study=="TN09"])
pca.TN09.rna_batch_removed <-
  calc_PCAs(vwts.TN09.rna_batch_removed, cpm=FALSE, log2_transform=FALSE)
scores_design_pca.TN09.rna_batch_removed <-
  merge(
    master.final[master.final$study=="TN09",],
    as.data.frame(
      pca.TN09.rna_batch_removed$x),
    by.x = "libid", by.y="row.names")

# now plot them
for (j in c("source", "rna_processing_core", "rna_batch"))
  plot_PCAs(
    scores_design_pca.TN09.source_and_rna_processing_core_removed,
    pvars.labs=pca.TN09.source_and_rna_processing_core_removed$pvars.labs, PCs=1:3,
    color_by_var=j,
    my_cols=colorblind_pal()(8),
    file_prefix="plots/T1D_placebos_TN09.source_and_rna_processing_core_removed",
    plotdims=c(9,6))

for (j in c("source", "rna_processing_core", "rna_batch"))
  plot_PCAs(
    scores_design_pca.TN09.rna_batch_removed,
    pvars.labs=pca.TN09.rna_batch_removed$pvars.labs, PCs=1:3,
    color_by_var=j,
    my_cols=colorblind_pal()(8),
    file_prefix="plots/T1D_placebos_TN09.rna_batch_removed",
    plotdims=c(9,6))

# use pcCors to check removal of batch effect
cor_all_clinvars_vs_pcaAll.counts.TN09.source_and_rna_processing_core_removed <-
  calc_PCcors(
    pca.TN09.source_and_rna_processing_core_removed,
    master.final[master.final$study=="TN09",])
try(plot_PCcor_heatmap(
  cor_all_clinvars_vs_pcaAll.counts.TN09.source_and_rna_processing_core_removed[1:10,],
  filename=
    "plots/T1D_plac.heatmap_correlations_clinvars_pcaAll.counts.TN09.source_and_rna_processing_core_removed.pdf",
  plotdims=c(12,6)))

cor_all_clinvars_vs_pcaAll.counts.TN09.rna_batch_removed <-
  calc_PCcors(
    pca.TN09.rna_batch_removed,
    master.final[master.final$study=="TN09",])
try(plot_PCcor_heatmap(
  cor_all_clinvars_vs_pcaAll.counts.TN09.rna_batch_removed[1:10,],
  filename=
    "plots/T1D_plac.heatmap_correlations_clinvars_pcaAll.counts.TN09.rna_batch_removed.pdf",
  plotdims=c(12,6)))

## for TN09, both work similarly well; use rna_batch


##### generate vwts.all.removed_study_rna_batch #####

DesignMat.removed_study_rna_batch <-
  model.matrix(
    ~ log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept +
      age_years + sex_by_rna,
    data=master.final)
vwts.all.removed_study_rna_batch <-
  removeBatchEffect(
    vwts.all[,match(master.final$libid, colnames(vwts.all))],
    batch=master.final$study_rna_batch,
    design=DesignMat.removed_study_rna_batch)


##### calculate median gene set expression for CXCR1.mod and CD19.mod #####

## read in Peter's gene sets (derived from Cate's immune gene atlas)
gene_sets.Linsley <-
  read_delim(
    "Supp_Table_1_Genes in immune molecular modules inc named genes.gmx",
    delim="\t") %>%
  as.list() %>% # convert to a list of character vectors of genes
  lapply(FUN=na.omit) # remove NAs

## calculate median expression after removing study_rna_batch, in log2(cpm+1) of CXCR1.mod, CD19.mod, and MPO.mod
# store in master.final
master.final$median_CXCR1.mod.removed_study_rna_batch <-
  gene_set_median_count(
    gene_sets.Linsley[["CXCR1.mod"]],
    vwts.all.removed_study_rna_batch[
      , match(master.final$libid,
              colnames(vwts.all.removed_study_rna_batch))])

master.final$median_CD19.mod.removed_study_rna_batch <-
  gene_set_median_count(
    gene_sets.Linsley[["CD19.mod"]],
    vwts.all.removed_study_rna_batch[
      , match(master.final$libid,
              colnames(vwts.all.removed_study_rna_batch))])

master.final$median_MPO.mod.removed_study_rna_batch <-
  gene_set_median_count(
    gene_sets.Linsley[["MPO.mod"]],
    vwts.all.removed_study_rna_batch[
      , match(master.final$libid,
              colnames(vwts.all.removed_study_rna_batch))])


##### compare gene set expression values to each other (are high neutrophil libs also high B cell libs) #####

# calculate correlation of median CXCR1.mod and CD19.mod expression
cor(
  master.final[
    master.final$rnaseq_baseline_visit,  # for baseline visits only, max of one per patient
    c("median_CXCR1.mod.removed_study_rna_batch", "median_CD19.mod.removed_study_rna_batch")],
  method="pearson")


##### merge values calculated from RNAseq back into master, for combination downstream with CBC data #####

master <-
  merge(master,
        master.final[
          ,c("libid", setdiff(colnames(master.final), colnames(master)))],
        by="libid", all.x=TRUE)


##### plot progressor status by age #####

# plot C-peptide rate vs. age, with logarithmic function fit
pdf("Fig_1C.pdf", w=9, h=6)
ggplot(
  data=master[!duplicated(master$participant_id),],
  mapping=aes(
    x=age_years,
    y=log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept)) +
  geom_point(size=2) +
  labs(x="Age at diagnosis", y="Rate of C-peptide change") +
  geom_smooth(method="lm", formula= y~I(log(x)), color="black")
dev.off()

# fit model to rate of C-peptide change vs. log age in years
lm.cpeptide_slope_vs_log_age_years <-
  lm(log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept ~ log(age_years),
     data = master %>%
       dplyr::filter(!duplicated(participant_id)))
summary(lm.cpeptide_slope_vs_log_age_years)

# test for heteroscedasticity of model fit to age in years
lmtest::bptest(
  lm(log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept ~ age_years,
     data = master %>%
       dplyr::filter(!duplicated(participant_id))))


##### plot AUC variables split out by age_groups or tertiles, or by age at draw #####

## plot log AUC by individual vs cpeptide_study_day
master.tmp <- master
baseline_auc.tmp <-
  master.tmp[master.tmp$cpeptide_baseline_visit,]

baseline_auc.tmp$baseline_auc_group <-
  cut_number(baseline_auc.tmp$log_auc2hr, n=3)
master.tmp <-
  merge(master.tmp,
        baseline_auc.tmp[,c("participant_id", "baseline_auc_group")],
        all.x=TRUE)


##### examine relationship between C-peptide slope and HLA risk score #####

## rate of T1D progression by HLA risk score (Winkler 2014)
pdf("Fig_1D.pdf", w=8, h=5)
ggplot(
  master %>%
    dplyr::filter(
      !duplicated(participant_id),
        !is.na(hla_risk_category_winkler_2014)),
  mapping=aes(
    x=hla_risk_category_winkler_2014,
    y=log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept)) +
  ggbeeswarm::geom_beeswarm(
    size=3, cex=2.5, shape=16) +
  stat_smooth(
    size=2, se=FALSE, method="lm", color="black", linetype="longdash", fullrange=TRUE) +
  scale_x_continuous(breaks=1:6, limits=c(0,7)) +
  coord_cartesian(xlim=c(0.6, 6.3)) +
  labs(x="HLA risk category", y="Rate of C-peptide change")
dev.off()

master %>%
  dplyr::filter(!duplicated(participant_id)) %>%
  lm(log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept ~ hla_risk_category_winkler_2014,
     data=.) %>%
  summary()

## rate of T1D progression by HLA risk score (Winkler 2014), age <= 15 years
pdf("Fig_S2.pdf", w=8, h=5)
ggplot(
  master %>%
    dplyr::filter(
      !duplicated(participant_id),
      !is.na(hla_risk_category_winkler_2014),
      age_years <= 15),
  mapping=aes(
    x=hla_risk_category_winkler_2014,
    y=log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept)) +
  ggbeeswarm::geom_beeswarm(
    size=3, cex=2.5, shape=16) +
  stat_smooth(
    size=2, se=FALSE, method="lm", color="black", linetype="longdash", fullrange=TRUE) +
  scale_x_continuous(breaks=1:6, limits=c(0,7)) +
  coord_cartesian(xlim=c(0.6, 6.3)) +
  labs(x="HLA risk category", y="Rate of C-peptide change")
dev.off()

master %>%
  dplyr::filter(!duplicated(participant_id), age_years <= 15) %>%
  lm(log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept ~ hla_risk_category_winkler_2014,
     data=.) %>%
  summary()


##### plot C-peptide levels over time, with individuals highlighted to illustrate variation #####

## plot log AUC by individual vs cpeptide_study_day
# with points and lines for two individuals
pdf("Fig_S1A.pdf", w=8, h=5.3)
ggplot(
  data=master[!is.na(master$log_auc2hr),],
  mapping=aes(x=cpeptide_study_day, y=log_auc2hr, group=participant_id)) +
  geom_polygon(
    inherit.aes=FALSE,
    data=data.frame(
      cpeptide_study_day=rep(range(with(master, diagnosis_date - day_0), na.rm=TRUE), each=2),
      log_auc2hr=
        (range(master$log_auc2hr, na.rm=TRUE) + c(-0.5,0.5))[c(1,2,2,1)]),
    mapping=aes(x=cpeptide_study_day, y=log_auc2hr),
    fill="gray80") +
  geom_line(size=0.3, color="gray70") +
  geom_point(
    data=master[master$participant_id=="TN02_9057024",],
    color="red", size=4) + # plot one fast progressor
  geom_smooth(
    data=master[
      with(master, participant_id=="TN02_9057024" & !prev_at_cpep_detect_thresh),],
    se=FALSE, method="lm", color="red") +
  geom_point(
    data=master[master$participant_id=="T1DAL_01710726",],
    color="blue", size=4) + # plot one slow progressor
  geom_smooth(
    data=master[
      with(master, participant_id=="T1DAL_01710726" & !prev_at_cpep_detect_thresh),],
    se=FALSE, method="lm", color="blue") +
  scale_x_continuous(breaks=seq(0,1000, by=250)) +
  scale_y_continuous(breaks=log(10^(-2:0)), labels=10^(-2:0)) +
  labs(x="Days in study", y="C-peptide 2hr AUC\n(nmol / L / min)") +
  geom_hline(yintercept=log(low_lim_detect), linetype="dashed") +
  coord_cartesian(
    xlim=xlim.cpeptide_AUC_figures,
    ylim=range(master$log_auc2hr, na.rm=TRUE))
dev.off()


##### summarize numbers of patients and samples #####

### total number of control-arm patients by study
# cpeptide_data.by_study
sapply(
  cpeptide_data.by_study,
  function(x) {
    x %>%
      dplyr::filter(treatment %in% "control") %>%
      dplyr::pull(participant_id) %>%
      unique() %>%
      length()})

# using cpeptide_data.merged
length(unique(cpeptide_data.merged$participant_id))
# 152 patients in the control cpeptide_data

# using patient_data.merged
length(unique(patient_data.merged$participant_id))
# 154 patients in the control-arm patient_data

patient_data.merged[
  patient_data.merged$participant_id %nin% unique(cpeptide_data.merged$participant_id),]
# 2 patients in AbATE who only had screening visits; leave them out


### C-peptide visits per patient
# using cpeptide_data.by_study
sapply(
  cpeptide_data.by_study,
  function(x) {
    x %>%
      dplyr::filter(treatment %in% "control") %>%
      dplyr::filter(!is.na(pep0)) %>%
      nrow()}) %>%
  sum() / length(unique(cpeptide_data.merged$participant_id))

# using cpeptide_data.merged
sapply(
  cpeptide_data.by_study,
  function(x) {
    x %>%
      dplyr::filter(treatment %in% "control") %>%
      dplyr::filter(!is.na(pep0)) %>%
      nrow()})

### total patients with RNA-seq data by study
# using rnaseq_annotation.by_study
sapply(
  rnaseq_annotation.by_study,
  function(x) {
    x %>%
      merge(patient_data.merged[,c("participant_id", "treatment")], all=FALSE) %>%
      dplyr::filter(treatment %in% "control") %>%
      dplyr::select(participant_id) %>%
      unique() %>%
      nrow()})


### number of RNA-seq visits by study

# using rnaseq_annotation.by_study
sapply(
  rnaseq_annotation.by_study,
  function(x) {
    x %>%
      merge(patient_data.merged[,c("participant_id", "treatment")], all=FALSE) %>%
      dplyr::filter(treatment %in% "control") %>%
      dplyr::select(rnaseq_visit_id) %>%
      unique() %>%
      nrow()})


##### summarize clinical values at baseline #####

### determine sex of patients by study in the C-peptide data
lapply(
  patient_data.by_study,
  function(x) {
    x %>%
      dplyr::filter(treatment %in% "control") %>%
      dplyr::filter(participant_id %in% cpeptide_data.merged$participant_id) %>%
      dplyr::select(participant_id, sex) %>%
      unique() %>%
      table() %>%
      colSums()})

### determine sex of patients by study in the C-peptide data
# using patient_data.by_study
lapply(
  patient_data.by_study,
  function(x) {
    x %>%
      dplyr::filter(treatment %in% "control") %>%
      dplyr::filter(participant_id %in% rnaseq_annotation.merged$participant_id) %>%
      dplyr::select(participant_id, sex) %>%
      unique() %>%
      table() %>%
      colSums()})


### summarize HbA1c at baseline by study

## for subjects with C-peptide data

# using cpeptide_data.merged for the ones I can (AbATE, T1DAL, TN02)
cpeptide_data.merged %>%
  dplyr::filter(cpeptide_baseline_visit & !is.na(hba1c)) %>%
  dplyr::select(study, participant_id, hba1c) %>%
  group_by(study) %>%
  summarise(mean=mean(hba1c), sd=sd(hba1c))

# using cpeptide_data.by_study for START (baselines: HbA1c at baseline visit, AUC at screening visit)
cpeptide_data.by_study[["START"]] %>%
  dplyr::filter(treatment == "control") %>%
  dplyr::filter(!is.na(baseline_hba1c)) %>%
  dplyr::select(participant_id, baseline_hba1c) %>%
  unique() %>%
  summarise(mean=mean(baseline_hba1c), sd=sd(baseline_hba1c))

# using cpeptide_data.by_study TN09 (baselines: HbA1c at baseline visit, AUC at screening visit)
cpeptide_data.by_study[["TN09"]] %>%
  dplyr::filter(treatment == "control") %>%
  dplyr::filter(cpeptide_visit_name=="baseline" & !is.na(hba1c)) %>%
  dplyr::select(participant_id, hba1c) %>%
  summarise(mean=mean(hba1c), sd=sd(hba1c))

# summarize for all subjects
(cpeptide_data.merged %>%
  dplyr::filter(cpeptide_baseline_visit & !is.na(hba1c)) %>%
  dplyr::select(participant_id, hba1c)) %>%
  rbind(
    cpeptide_data.by_study[["START"]] %>%
      dplyr::filter(treatment == "control") %>%
      dplyr::filter(!is.na(baseline_hba1c)) %>%
      dplyr::select(participant_id, baseline_hba1c) %>%
      unique() %>%
      dplyr::rename(hba1c=baseline_hba1c)) %>%
  rbind(
    cpeptide_data.by_study[["TN09"]] %>%
      dplyr::filter(treatment == "control") %>%
      dplyr::filter(cpeptide_visit_name=="baseline" & !is.na(hba1c)) %>%
      dplyr::select(participant_id, hba1c)) %>%
  summarise(mean=mean(hba1c), sd=sd(hba1c))


## for subjects with RNA-seq data

# using cpeptide_data.merged for the ones I can (AbATE, T1DAL, TN02)
cpeptide_data.merged %>%
  dplyr::filter(participant_id %in% rnaseq_annotation.merged$participant_id) %>%
  dplyr::filter(cpeptide_baseline_visit & !is.na(hba1c)) %>%
  dplyr::select(study, participant_id, hba1c) %>%
  group_by(study) %>%
  summarise(mean=mean(hba1c), sd=sd(hba1c))

# using cpeptide_data.by_study for START (baselines: HbA1c at baseline visit, AUC at screening visit)
cpeptide_data.by_study[["START"]] %>%
  dplyr::filter(participant_id %in% rnaseq_annotation.merged$participant_id) %>%
  dplyr::filter(treatment == "control") %>%
  dplyr::filter(!is.na(baseline_hba1c)) %>%
  dplyr::select(participant_id, baseline_hba1c) %>%
  unique() %>%
  summarise(mean=mean(baseline_hba1c), sd=sd(baseline_hba1c))

# using cpeptide_data.by_study TN09 (baselines: HbA1c at baseline visit, AUC at screening visit)
cpeptide_data.by_study[["TN09"]] %>%
  dplyr::filter(participant_id %in% rnaseq_annotation.merged$participant_id) %>%
  dplyr::filter(treatment == "control") %>%
  dplyr::filter(cpeptide_visit_name=="baseline" & !is.na(hba1c)) %>%
  dplyr::select(participant_id, hba1c) %>%
  summarise(mean=mean(hba1c), sd=sd(hba1c))

# summarize for all subjects
(cpeptide_data.merged %>%
    dplyr::filter(cpeptide_baseline_visit & !is.na(hba1c)) %>%
    dplyr::select(participant_id, hba1c)) %>%
  rbind(
    cpeptide_data.by_study[["START"]] %>%
      dplyr::filter(treatment == "control") %>%
      dplyr::filter(!is.na(baseline_hba1c)) %>%
      dplyr::select(participant_id, baseline_hba1c) %>%
      unique() %>%
      dplyr::rename(hba1c=baseline_hba1c)) %>%
  rbind(
    cpeptide_data.by_study[["TN09"]] %>%
      dplyr::filter(treatment == "control") %>%
      dplyr::filter(cpeptide_visit_name=="baseline" & !is.na(hba1c)) %>%
      dplyr::select(participant_id, hba1c)) %>%
  dplyr::filter(participant_id %in% rnaseq_annotation.merged$participant_id) %>%
  summarise(mean=mean(hba1c), sd=sd(hba1c))


### summarize age by study

## for subjects with C-peptide data

# using patient_data.merged
patient_data.merged %>%
  dplyr::filter(participant_id %in% cpeptide_data.merged$participant_id) %>%
  dplyr::select(study, age_years) %>%
  group_by(study) %>%
  summarise(median=median(age_years), min=min(age_years), max=max(age_years))

# using patient_data.merged for entire study
patient_data.merged %>%
  dplyr::filter(participant_id %in% cpeptide_data.merged$participant_id) %>%
  dplyr::select(age_years) %>%
  summarise(median=median(age_years), min=min(age_years), max=max(age_years))


## for subjects with RNA-seq data

# using patient_data.merged
patient_data.merged %>%
  dplyr::filter(participant_id %in% rnaseq_annotation.merged$participant_id) %>%
  dplyr::select(study, age_years) %>%
  group_by(study) %>%
  summarise(median=median(age_years), min=min(age_years), max=max(age_years))

# using patient_data.merged for entire study
patient_data.merged %>%
  dplyr::filter(participant_id %in% rnaseq_annotation.merged$participant_id) %>%
  dplyr::select(age_years) %>%
  summarise(median=median(age_years), min=min(age_years), max=max(age_years))


### summarize baseline C-peptide 2-hour AUC by study (in nmol/L)

## for subjects with C-peptide data

# using cpeptide_data.merged for the ones I can
cpeptide_data.merged %>%
  dplyr::filter(cpeptide_baseline_visit & !is.na(auc2hr)) %>%
  dplyr::select(study, participant_id, auc2hr) %>%
  group_by(study) %>%
  summarise(mean=mean(auc2hr), sd=sd(auc2hr))

# using cpeptide_data.merged for entire study
cpeptide_data.merged %>%
  dplyr::filter(cpeptide_baseline_visit & !is.na(auc2hr)) %>%
  dplyr::select(study, participant_id, auc2hr) %>%
  summarise(mean=mean(auc2hr), sd=sd(auc2hr))


## for subjects with RNA-seq data

# using cpeptide_data.merged for the ones I can
cpeptide_data.merged %>%
  dplyr::filter(participant_id %in% rnaseq_annotation.merged$participant_id) %>%
  dplyr::filter(cpeptide_baseline_visit & !is.na(auc2hr)) %>%
  dplyr::select(study, participant_id, auc2hr) %>%
  group_by(study) %>%
  summarise(mean=mean(auc2hr), sd=sd(auc2hr))

# using patient_data.merged for entire study
cpeptide_data.merged %>%
  dplyr::filter(participant_id %in% rnaseq_annotation.merged$participant_id) %>%
  dplyr::filter(cpeptide_baseline_visit & !is.na(auc2hr)) %>%
  dplyr::select(study, participant_id, auc2hr) %>%
  summarise(mean=mean(auc2hr), sd=sd(auc2hr))


##### predict 2-year C-peptide AUC from baseline and 6-month C-peptide AUC and/or age #####

## models to use
# baseline C-peptide
# baseline C-peptide + age
# baseline C-peptide + 6-month C-peptide
# baseline C-peptide + 12-month C-peptide

cpeptide_visit_name.6mo.tmp <-
  c("month 6", "26 weeks", "week 24", "day 196")
cpeptide_visit_name.12mo.tmp <-
  c("month 12", "52 weeks", "week 52", "day 364")
cpeptide_visit_name.24mo.tmp <-
  c("month 24", "104 weeks", "week 104", "day 728")

cpeptide_data.merged$log_auc2hr <-
  log(cpeptide_data.merged$auc2hr)

cpeptide_data.merged.for_model <-
  unique(
    cpeptide_data.merged[
      cpeptide_data.merged$cpeptide_baseline_visit,
      c("participant_id", "age_years", "log_auc2hr")]) %>%
  dplyr::rename(log_auc2hr_baseline=log_auc2hr)
cpeptide_data.merged.for_model$log_auc2hr_6mo <-
  cpeptide_data.merged$log_auc2hr[
    cpeptide_data.merged$cpeptide_visit_name %in%
      cpeptide_visit_name.6mo.tmp][
        match(cpeptide_data.merged.for_model$participant_id,
              cpeptide_data.merged$participant_id[
                cpeptide_data.merged$cpeptide_visit_name %in%
                  cpeptide_visit_name.6mo.tmp])]
cpeptide_data.merged.for_model$log_auc2hr_12mo <-
  cpeptide_data.merged$log_auc2hr[
    cpeptide_data.merged$cpeptide_visit_name %in%
      cpeptide_visit_name.12mo.tmp][
        match(cpeptide_data.merged.for_model$participant_id,
              cpeptide_data.merged$participant_id[
                cpeptide_data.merged$cpeptide_visit_name %in%
                  cpeptide_visit_name.12mo.tmp])]
cpeptide_data.merged.for_model$log_auc2hr_24mo <-
  cpeptide_data.merged$log_auc2hr[
    cpeptide_data.merged$cpeptide_visit_name %in%
      cpeptide_visit_name.24mo.tmp][
        match(cpeptide_data.merged.for_model$participant_id,
              cpeptide_data.merged$participant_id[
                cpeptide_data.merged$cpeptide_visit_name %in%
                  cpeptide_visit_name.24mo.tmp])]

cpeptide_data.merged.for_model <-
  cpeptide_data.merged.for_model[complete.cases(cpeptide_data.merged.for_model),]
# 109 complete cases

formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years <-
  list(formula(log_auc2hr_24mo ~ log_auc2hr_baseline),
       formula(log_auc2hr_24mo ~ log_auc2hr_baseline + age_years),
       formula(log_auc2hr_24mo ~ log_auc2hr_baseline + log_auc2hr_6mo),
       formula(log_auc2hr_24mo ~ log_auc2hr_baseline + log_auc2hr_12mo)) %>%
  setNames(
    c("log_auc2hr_baseline",
      "log_auc2hr_baseline_age_years", 
      "log_auc2hr_baseline_6mo",
      "log_auc2hr_baseline_12mo"))
tobit.formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years <-
  lapply(
    formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years,
    function(x) tobit(x, data=cpeptide_data.merged.for_model, left=log(low_lim_detect)))


#### leave-one-out cross-validation of the models

## need to do manually because there aren't built-ins for the tobit model

## store the predicted values (for plotting)
pred_LOOCV.tobit.formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years.test <-
  list()

## store the PRESS statistic values
press_LOOCV.tobit.formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years <-
  replicate(
    n=length(formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years),
    as.numeric(rep(NA, 10)), simplify=FALSE) %>%
  setNames(names(formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years))

## iterate over individuals
n <- nrow(cpeptide_data.merged.for_model)
for (i in (1:n)) {
  
  ## select training and test data
  cpeptide_data.merged.for_model.training.tmp <-
    cpeptide_data.merged.for_model[-i,]
  cpeptide_data.merged.for_model.test.tmp <-
    cpeptide_data.merged.for_model[i,]
  
  ## fit model on all but one observations
  tobit_LOOCV.formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years.training.tmp <-
    lapply(
      formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years,
      function(x)
        tobit(x, data=cpeptide_data.merged.for_model.training.tmp,
              left=log(low_lim_detect)))
  
  ## predict for remaining observations
  pred_LOOCV.tobit.formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years.test[[i]] <-
    lapply(
      tobit_LOOCV.formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years.training.tmp,
      function(x) {
        predict(
          x,
          newdata=
            cpeptide_data.merged.for_model.test.tmp) %>%
          sapply(function(y) max(y, log(low_lim_detect)))})
  
  ## calculate PRESS statistic (predicted residual sum of squares)
  for (j in names(formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years)) {
    press_LOOCV.tobit.formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years[[j]][i] <-
      sum(
        (cpeptide_data.merged.for_model.test.tmp$log_auc2hr_24mo -
           pred_LOOCV.tobit.formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years.test[[i]][[j]])^2)}
}  

## calculate sum PRESS statistic across all observations
lapply(press_LOOCV.tobit.formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years, sum)
press_LOOCV.tobit.formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years.df <-
  data.frame(
    model=
      rep(names(press_LOOCV.tobit.formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years),
          each=n),
    press=
      unlist(press_LOOCV.tobit.formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years))

### make plots of predicted vs. observed, using prediction from each leave-one-out
pred_obs_LOOCV.tobit.formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years <-
  list()
for (i in names(formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years)) {
  pred_obs_LOOCV.tobit.formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years[[i]] <-
    data.frame(
      predicted=
        unlist(
          lapply(pred_LOOCV.tobit.formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years.test,
                 function(x) x[[i]])),
      observed=
        cpeptide_data.merged.for_model$log_auc2hr_24mo)
}

### output plots of predicted vs. observed for all samples from all CV folds
lims.tmp <-
  pred_obs_LOOCV.tobit.formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years %>%
  lapply(function(x) x[,c("predicted", "observed")]) %>%
  unlist() %>%
  range(na.rm=TRUE)

for (model.tmp in names(formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years)) {
  pdf(paste0("Fig_1B_", model.tmp, ".LOO_CV.pdf"), w=4.3, h=4)
  print(
    ggplot(
      pred_obs_LOOCV.tobit.formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years[[model.tmp]],
      mapping=aes(
        x=predicted,
        y=observed)) +
      geom_point(size=2.5, alpha=0.8) +
      geom_abline(size=1, slope=1, intercept=0, linetype="dashed") +
      scale_x_continuous(
        breaks=log(10^(-2:0)), labels=10^(-2:0),
        limits=lims.tmp) +
      scale_y_continuous(
        breaks=log(10^(-2:0)), labels=10^(-2:0),
        limits=lims.tmp) +
      labs(x="predicted 24-month AUC", y="observed 24-month AUC", title=model.tmp) +
      annotate(
        geom="text",
        x=0.5, y=-3.5,
        size=4,
        label=
          paste(
            "pred R^2 =",
            round(cor(pred_obs_LOOCV.tobit.formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years[[model.tmp]][,"predicted"],
                      pred_obs_LOOCV.tobit.formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years[[model.tmp]][,"observed"],
                      use="pairwise")^2,
                  2))) +
      annotate(
        geom="text",
        x=0.5, y=-4,
        size=4,
        label=
          paste(
            "PRESS =",
            round(sum(press_LOOCV.tobit.formula_list.log_auc2hr_24mo.baseline_6mo_12mo_age_years[[model.tmp]]),
                  1))))
  dev.off()
}

rm_tmp(ask=FALSE)


##### save objects for downstream use #####

save(file="T1D_placebos_data_1_for_downstream_analyses.RData",
     list=c(
       "patient_data.merged", "patient_data.by_study",
       "study_schedules",
       "cbc.merged",
       "cpeptide_data.merged",
       "rnaseq_annotation.merged",
       "counts.final", "vwts.all", "vwts.all.removed_study_rna_batch",
       "master", "master.final",
       "gene_sets.Linsley"))
