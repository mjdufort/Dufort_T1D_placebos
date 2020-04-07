##### scripts for analysis of data for Dufort et al. 2019. Cell type–specific immune phenotypes predict loss of insulin secretion in new-onset type 1 diabetes. (DOI: 10.1172/jci.insight.125556)
### this file includes scripts for analyses of the data from the TN-05 rituximab trial

##### set up environment: load packages #####

## load general packages
library(tidyverse)
theme_set(
  theme_bw(20) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour="black", fill=NA, size=1),
          axis.text=element_text(colour="black"),
          axis.ticks=element_line(colour="black")))
update_geom_defaults("point", list(shape=16))
library(ggthemes)

## load analysis-specific packages
library(limma)
library(edgeR)
library(lmerTest)
library(annotables)
library(datamart)
library(ordinal)

# load useful functions
library(RNAseQC)
library(countSubsetNorm)
library(miscHelpers)
library(geneSetTools)
library(limmaTools)

# load custom packages (available at github.com/benaroyaresearch)
if (!require(RNAseQC)) remotes::install_github("benaroyaresearch/RNAseQC"); library(RNAseQC)
if (!require(countSubsetNorm)) remotes::install_github("benaroyaresearch/countSubsetNorm"); library(countSubsetNorm)
if (!require(miscHelpers)) remotes::install_github("benaroyaresearch/miscHelpers"); library(miscHelpers)
if (!require(geneSetTools)) remotes::install_github("benaroyaresearch/geneSetTools"); library(geneSetTools)
if (!require(limmaTools)) remotes::install_github("benaroyaresearch/limmaTools"); library(limmaTools)


##### set up color palette #####

color.placebo <- ggthemes::colorblind_pal()(8)[7]
color.ritux <- ggthemes::colorblind_pal()(8)[6]


##### load/save data from previous scripts #####

load("T1D_placebos_data_1_for_downstream_analyses.RData")
load("T1D_placebos_data_2_for_downstream_analyses.RData")
load("T1D_placebos_data_3_for_downstream_analyses.RData")
load("T1D_placebos_data_4_for_downstream_analyses.RData")


##### load patient and clinical data #####

## patient data (including age, sex, HLA)
patient_data_TN05 <-
  read.csv("TN05_patient_data.csv", stringsAsFactors = FALSE) %>%
  mutate(
    treatment = factor(treatment, levels=c("placebo", "rituximab")),
    age_split2 = factor(age_split2, levels=c("[8.34,16.1]", "(16.1,40.5]")))

## CBC data
cbc_data_TN05 <- read.csv("TN05_cbc_data_processed_merged.csv", stringsAsFactors = FALSE)
# CBC data have already had absolute numbers log-transformed and percentages arcsin-square-root transformed

## C-peptide data
cpeptide_data_TN05 <-
  read.csv("TN05_cpeptide_data_processed_merged.csv", stringsAsFactors = FALSE) %>%
  left_join(
    patient_data_TN05 %>% dplyr::select(study, participant_id, treatment, age_split2))


##### calculate C-peptide derivative values #####

## determine baseline values, and calculate log auc
cols.tmp <-
  c("baseline_auc2hr_nmol_l_min", "baseline_log_auc2hr_nmol_l_min")
cpeptide_data_TN05[, cols.tmp] <- as.numeric(NA)

for (i in unique(cpeptide_data_TN05$participant_id)) {
  rows.tmp <- which(cpeptide_data_TN05$participant_id==i)
  data.tmp <- cpeptide_data_TN05[rows.tmp,]
  baseline_row.tmp <- which(data.tmp$cpeptide_visit_name %in% c("baseline", "0"))

  if (length(baseline_row.tmp)==0)
    baseline_row.tmp <- which(data.tmp$cpeptide_visit_name %in% c("screening", "-1"))

  if (length(baseline_row.tmp)==0) next

  cpeptide_data_TN05[rows.tmp, cols.tmp] <-
    data.tmp[baseline_row.tmp, str_replace(cols.tmp, "baseline_", "")]
}

rm_tmp(ask=FALSE)


##### determine visits with C-peptide previously at detection threshold, to exclude them #####

cpeptide_data_TN05$prev_at_cpep_detect_thresh <- as.logical(NA)
cpeptide_data_TN05 <-
  cpeptide_data_TN05 %>%
  arrange(participant_id, cpeptide_study_day)

for (i in unique(cpeptide_data_TN05$participant_id)) { # iterate over patient
  rows.tmp <-
    which(
      with(cpeptide_data_TN05,
           participant_id==i & !is.na(auc2hr_nmol_l_min))) # determine all rows for current patient
  cpeptide_data_TN05$prev_at_cpep_detect_thresh[rows.tmp[1]] <- FALSE # set first visit to FALSE
  if (length(rows.tmp) > 1)
    for (j in (2:length(rows.tmp))) {
      cpeptide_data_TN05$prev_at_cpep_detect_thresh[rows.tmp[j]] <-
        (cpeptide_data_TN05$auc2hr_nmol_l_min[rows.tmp[j-1]] < low_lim_detect * 1.01) |  # previous visit is at detection limit
        cpeptide_data_TN05$prev_at_cpep_detect_thresh[rows.tmp[j-1]]  # previous visit had prev_at_cpep_detect_thresh
    }
}

rm_tmp(ask=FALSE)


##### fit curves to c-peptide AUC by cpeptide_study_day #####

# scale cpeptide_study_day so that it works better in model fitting
cpeptide_study_day_scaled.tmp <-
  scale(cpeptide_data_TN05$cpeptide_study_day, center=FALSE,
        scale=sd(cpeptide_data_TN05$cpeptide_study_day, na.rm=TRUE))
cpeptide_data_TN05$cpeptide_study_day_scaled <-
  as.vector(cpeptide_study_day_scaled.tmp)
days.sd <- attr(cpeptide_study_day_scaled.tmp, "scaled:scale")

## fit a mixed-effects model with a term to allow different fixed-effect slopes by treatment
lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random <-
  lmer(
    log_auc2hr_nmol_l_min ~
      (1|participant_id) +
      cpeptide_study_day_scaled +
      cpeptide_study_day_scaled:treatment +
      (cpeptide_study_day_scaled-1|participant_id),
    data=cpeptide_data_TN05[!cpeptide_data_TN05$prev_at_cpep_detect_thresh,])
summary(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)


##### extract slopes from model fit to log_auc2hr by cpeptide_study_day #####

cpeptide_data_TN05$log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept <-
  as.numeric(NA)

for (i in 1:nrow(cpeptide_data_TN05)) {
  # extract slope for linear random with intercept, with separate fixed effects for placebo and rituximab
  if (sum(!is.na(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id[
    match(
      cpeptide_data_TN05$participant_id[i],
      rownames(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id)),])) == 3)
    cpeptide_data_TN05$log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept[i] <-
      (coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id$cpeptide_study_day_scaled[
        match(cpeptide_data_TN05$participant_id[i],
              rownames(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id))] +
         coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id$'cpeptide_study_day_scaled:treatmentrituximab'[
           match(cpeptide_data_TN05$participant_id[i],
                 rownames(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id))] *
         (cpeptide_data_TN05$treatment[i]=="rituximab")) /
      days.sd * 365.25
}

rm_tmp(ask=FALSE)


##### plot and model rate of C-peptide change vs. age split and treatment, with interaction #####

## plot C-peptide rate
pdf("Fig_S12.pdf", w=11, h=6)
ggplot(
  cpeptide_data_TN05[!duplicated(cpeptide_data_TN05$participant_id),],
  mapping=aes(
    x=age_split2,
    y=log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept,
    fill=treatment)) +
  geom_boxplot(outlier.color=NA, position=position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(jitter.width=0.1, dodge.width=0.8)) +
  scale_fill_manual(
    values=c("rituximab"=make_transparent(color.ritux, alpha=150), "placebo"=color.placebo)) +
  labs(x="Age Group", y="Rate of C-peptide change") +
  scale_x_discrete(
    labels=c("[8.34,16.1]"="Younger half\n(8.3 to 16.1 years)",
             "(16.1,40.5]"="Older half\n(>16.1 years)"))
dev.off()

## fit model to C-peptide slope vs. age split and treatment, with interaction
lm.cpeptide_slope_vs_age_split2.by_treatment <-
  lm(data=cpeptide_data_TN05[!duplicated(cpeptide_data_TN05$participant_id),],
     formula=
       log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept~
       age_split2 + treatment + age_split2:treatment)
summary(lm.cpeptide_slope_vs_age_split2.by_treatment)

contrasts.tmp <-
  matrix(c(0,0,1,1), # rituximab effect in older subjects
         nrow=1, byrow=TRUE)
multcomp::glht(
  lm.cpeptide_slope_vs_age_split2.by_treatment,
  linfct=contrasts.tmp, alternative="two.sided") %>%
  summary()


##### pull in RNAseq data, and match with other data #####

## read in quality-controled, filtered, normalized counts
# split into two files to facilitate storage on github
counts_TN05 <-
  cbind(
    read.table(
      "TN05_counts_normalized_1.txt",
      sep="\t", header=T),
    read.table(
      "TN05_counts_normalized_2.txt",
      sep="\t", header=T)) %>%
  # fix Excel-mangled gene names
  magrittr::set_rownames(
    rownames(.) %>%
      str_replace("^(?=[0-9]+\\-Mar)", "MARCH") %>%
      str_replace("\\-Mar$", "") %>%
      str_replace("^(?=[0-9]+\\-Sep)", "SEPT") %>%
      str_replace("\\-Sep$", ""))


## read in sample annotation
rnaseq_annotation_TN05 <- read.csv("TN05_rnaseq_annotation.csv", stringsAsFactors = FALSE)

# merge patient_data_TN05 into rnaseq_annotation_TN05 data
rnaseq_annotation_TN05 <-
  merge(rnaseq_annotation_TN05,
        patient_data_TN05[
          c("participant_id",
            setdiff(colnames(patient_data_TN05), colnames(rnaseq_annotation_TN05)))],
        by="participant_id", all.x=T)

# generate visit_id
rnaseq_annotation_TN05$rnaseq_visit_id <-
  with(rnaseq_annotation_TN05,
       paste(participant_id, rnaseq_visit_number, sep="_"))
rnaseq_annotation_TN05$cpeptide_visit_id <-
  with(rnaseq_annotation_TN05,
       paste(participant_id, cpeptide_visit_number, sep="_"))

# merge cpeptide_data_TN05 into rnaseq_annotation_TN05
rnaseq_annotation_TN05 <-
  merge(rnaseq_annotation_TN05,
        cpeptide_data_TN05[
          c("cpeptide_visit_id",
            setdiff(colnames(cpeptide_data_TN05), colnames(rnaseq_annotation_TN05)))],
        by="cpeptide_visit_id", all.x=T)

rm_tmp(ask=FALSE)


##### filter rnaseq_annotation_TN05 and counts_TN05 to include only shared libraries, and put them in the same order #####

rnaseq_annotation_TN05 <-
  rnaseq_annotation_TN05[rnaseq_annotation_TN05$libid %in% colnames(counts_TN05),] %>%
  dplyr::arrange(participant_id, rnaseq_visit_number)
counts_TN05 <- counts_TN05[,match(rnaseq_annotation_TN05$libid, colnames(counts_TN05))]
# dim(rnaseq_annotation_TN05) # 197 libraires


##### some basic checks for problematic data #####

# calculate ratio of reads mapping to X and to Y chromosome, for each library
logXY.tmp <- logXYratio(counts_TN05, lib_cols=1:ncol(counts_TN05), gene_ID="symbol")
hist(logXY.tmp, breaks=15, main="log-ratio of X reads to Y reads")
cbind(rnaseq_annotation_TN05$sex,
      logXY.tmp[match(rnaseq_annotation_TN05$libid, names(logXY.tmp))])[
        order(logXY.tmp[match(rnaseq_annotation_TN05$libid, names(logXY.tmp))]),]
sort(logXY.tmp)
logXY_threshold <- 4

rnaseq_annotation_TN05$sex_by_rna <-
  ifelse(logXY.tmp[rnaseq_annotation_TN05$libid] > logXY_threshold, "F", "M")

# check that inferred sex is same for all libraries for each patient
table(rnaseq_annotation_TN05[,c("participant_id", "sex_by_rna")])
which(rowSums(table(rnaseq_annotation_TN05[,c("participant_id", "sex_by_rna")])==0)!=1)
# looks good!

# check that inferred sex matches sex from database
stopifnot(all(rnaseq_annotation_TN05$sex == rnaseq_annotation_TN05$sex_by_rna))
# looks good!

# drop problematic libraries, and create quality-controlled rnaseq_annotation_TN05 and counts_TN05 objects
# except no problematic libraries to drop here
rnaseq_annotation_TN05.qc <- rnaseq_annotation_TN05
counts_TN05.qc <- counts_TN05[,match(rnaseq_annotation_TN05.qc$libid, colnames(counts_TN05))]

rm_tmp(ask=FALSE)


##### generate simply-named combined objects with all libraries passing QC and outlier cuts #####

counts_TN05.final <- counts_TN05.qc
master_TN05.final <- rnaseq_annotation_TN05.qc

# split subjects by age above/below median
age_split.tmp <- unique(master_TN05.final[,c("participant_id", "age_years")])
age_split.tmp$age_split_rnaseq_only <- cut_number(age_split.tmp$age_years, n=2)
master_TN05.final$age_split_rnaseq_only <-
  age_split.tmp$age_split_rnaseq_only[
    match(master_TN05.final$participant_id, age_split.tmp$participant_id)]
master_TN05.final[,c("age_split_rnaseq_only", "age_years")] # looks good!

### merge CBC data into master_TN05.final
master_TN05.final <-
  base::merge(
    master_TN05.final,
    cbc_data_TN05[,setdiff(colnames(cbc_data_TN05), colnames(master_TN05.final))],
    by.x="rnaseq_visit_id",
    by.y="cbc_visit_id",
    all.x=TRUE) %>%
  arrange(rnaseq_visit_id)


##### calculate median expression of CD19.mod #####

vwts_all_TN05 <-
  calc_norm_counts(
    counts=counts_TN05.final, design=master_TN05.final, libID_col="libid",
    min_cpm=1, min_libs_perc=0.15,
    normalize=FALSE,
    return_DGEcounts=TRUE) %>%
  voom()

## store median gene set expression in master_TN05.final
master_TN05.final$median_CD19.mod <-
  gene_set_median_count(
    gene_sets.Linsley[["CD19.mod"]], vwts_all_TN05,
    remove_low_count_genes=TRUE)[
      match(master_TN05.final$libid, colnames(vwts_all_TN05$E))]

## determine baseline median scores for CD19.mod
master_TN05.final[, "baseline_median_CD19.mod"] <- as.numeric(NA)
for (i in unique(master_TN05.final$participant_id)) {
  master.tmp <-
    master_TN05.final[
      with(master_TN05.final,
           (participant_id==i) & (rnaseq_visit_number==2)),]
  if (nrow(master.tmp)==1) {
    master_TN05.final[
      master_TN05.final$participant_id==i,
      "baseline_median_CD19.mod"] <-
      master.tmp[, "median_CD19.mod"]
  } else cat(nrow(master.tmp), "baseline visits for", i, "\n")
}

##### plot and model rate of C-peptide change vs. median gene set expression, by treatment group #####

## plot it
pdf("Fig_5A.pdf", w=9, h=6)
print(
  ggplot(
    master_TN05.final[
      with(master_TN05.final,
           rnaseq_visit_number==2),],
    mapping=aes(
      x=median_CD19.mod,
      y=log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept,
      color=treatment,
      group=treatment)) +
    geom_smooth(method="lm", se=FALSE) +
    geom_point(size=3) +
    scale_color_manual(values=c("rituximab"=color.ritux, "placebo"=color.placebo)) +
    labs(x="B cell genes at baseline\n(median CD19.mod expression)", y="Rate of C-peptide change",
         title=NULL) +
    coord_cartesian(
      xlim=range(master_TN05.final$median_CD19.mod[master_TN05.final$rnaseq_visit_number==2], na.rm=T),
      ylim=range(
        master_TN05.final$log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept[
          master_TN05.final$rnaseq_visit_number==2], na.rm=T)))
dev.off()

## fit model to AUC slope by median_CD19.mod, by treatment
lm.cpeptide_slope_vs_baseline_median_CD19.mod.by_treatment <-
  lm(data=master_TN05.final[master_TN05.final$rnaseq_visit_number==2,],
     formula=
       log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept ~
       median_CD19.mod * treatment)
summary(lm.cpeptide_slope_vs_baseline_median_CD19.mod.by_treatment)

# determine slope for rituxiamb group
contrasts.tmp <-
  matrix(c(0,1,0,1),
         nrow=1, byrow=TRUE)
multcomp::glht(
  lm.cpeptide_slope_vs_baseline_median_CD19.mod.by_treatment,
  linfct=contrasts.tmp, alternative="two.sided") %>%
  summary()

rm_tmp(ask=FALSE)


##### plot and model rate of C-peptide change vs. B cell gene set expression, by age split #####

# plot rate of C-peptide change vs. B cell gene set expression, by treatment, for each age split
pdf("Fig_5C.pdf",
    w=9, h=6, useDingbats=FALSE)
translate_age.tmp <-
  c("Younger half", "Older half") %>%
  setNames(c("[12.1,19.1]", "(19.1,40.5]"))
for (i in unique(master_TN05.final$age_split_rnaseq_only)) {
  name.tmp <- i %>%
    str_replace_all("\\[|\\]|\\(|\\)", "") %>%
    str_replace_all(",Inf", " up") %>%
    str_replace_all(",", " to ")
  print(
    ggplot(
      master_TN05.final[
        with(master_TN05.final,
             rnaseq_visit_number==2 &
               age_split_rnaseq_only==i),],
      mapping=aes(
        x=median_CD19.mod,
        y=log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept,
        color=treatment,
        group=treatment)) +
      geom_smooth(method="lm", se=FALSE) +
      geom_point(size=3) +
      scale_color_manual(values=c("rituximab"=color.ritux, "placebo"=color.placebo)) +
      labs(x="B cell gene expression\n(median CD19.mod)",
           y="Rate of C-peptide change",
           title=paste0(translate_age.tmp[i], " (", name.tmp, " years)")) +
      coord_cartesian(
        xlim=range(master_TN05.final$median_CD19.mod[master_TN05.final$rnaseq_visit_number==2], na.rm=T),
        ylim=range(
          master_TN05.final$log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept[
            master_TN05.final$rnaseq_visit_number==2], na.rm=T)))
}
dev.off()

## fit model to cpeptide_slope, with median_CD19.mod, age_split_rnaseq_only, treatment, and interaction
lm.cpeptide_slope_vs_baseline_median_CD19.mod_age_split_rnaseq_only_treatment.with_interaction <-
  lm(data=master_TN05.final[master_TN05.final$rnaseq_visit_number==2,],
     formula=
       log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept~
       median_CD19.mod*age_split_rnaseq_only*treatment)
summary(lm.cpeptide_slope_vs_baseline_median_CD19.mod_age_split_rnaseq_only_treatment.with_interaction)
# slope for median_CD19.mod in placebo young patients is less than in rituximab young patients (p=0.0006; median_CD19.mod:treatmentrituximab)

# one-sided significance tests
contrasts.tmp <-
  matrix(c(0,0,0,0,0,1,0,0, # for placebo vs. rituximab in age_split_rnaseq_only_lo
           0,0,0,0,-1,1,0,1), # for placebo vs. rituximab in age_split_rnaseq_only_hi
         nrow=2, byrow=TRUE)
multcomp::glht(
  lm.cpeptide_slope_vs_baseline_median_CD19.mod_age_split_rnaseq_only_treatment.with_interaction,
  linfct=contrasts.tmp, alternative="greater") %>%
  summary()

contrasts.tmp <-
  matrix(c(0,0,0,0,1,0,0,-1), # for the difference between placebo/rituximab slope differences in the older and younger subjects
         nrow=1, byrow=TRUE)
multcomp::glht(
  lm.cpeptide_slope_vs_baseline_median_CD19.mod_age_split_rnaseq_only_treatment.with_interaction,
  linfct=contrasts.tmp, alternative="greater") %>%
  summary()

rm_tmp(ask=FALSE)


##### load flow data, and merge with other datasets #####

flow_data_TN05 <- read.csv("TN05_flow_data_processed_merged.csv", stringsAsFactors = FALSE)

## merge patient_data_TN05 into flow_data_TN05
flow_data_TN05 <-
  merge(flow_data_TN05, patient_data_TN05,
        all.x=TRUE)


##### merge flow_data with master RNAseq annotation object #####

intersect(colnames(flow_data_TN05), colnames(master_TN05.final))
master_flow_data_merged_TN05 <-
  merge(master_TN05.final,
        flow_data_TN05[
          ,c("study", "participant_id", "treatment",
             setdiff(colnames(flow_data_TN05), colnames(master_TN05.final)))],
        by.x=c("study", "participant_id", "treatment", "rnaseq_visit_number"),
        by.y=c("study", "participant_id", "treatment", "study_visit_number"),
        all=TRUE) %>%
  mutate(treatment = factor(treatment, levels=c("placebo", "rituximab")))

# check for all-NA rows
ncol(master_flow_data_merged_TN05)
rowSums(!is.na(master_flow_data_merged_TN05)) %>%
  min()
# none; one with 9 non-NA columns

# propagate patient-specific values (e.g. progressor status and C-peptide slope) to all datapoints for each patient
cols.to_propagate.tmp <-
  c("participant_id", "treatment", "age_years",
    "race", "ethnicity", "sex", "age_split2", "age_split_rnaseq_only",
    "baseline_auc2hr_nmol_l_min", "baseline_log_auc2hr_nmol_l_min",
    "log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept")

colSums(is.na(master_flow_data_merged_TN05))
colSums(is.na(master_flow_data_merged_TN05[,cols.to_propagate.tmp])) %>%
  sum()

for (i in cols.to_propagate.tmp) {
  rows.tmp <- which(is.na(master_flow_data_merged_TN05[,i]))
  cat("Rows to fill for variable", i, ": ", rows.tmp, "\n")
  for (j in rows.tmp) {
    # first try to get data from within data frame
    values.tmp <-
      master_flow_data_merged_TN05[
        match(master_flow_data_merged_TN05$participant_id[j], master_flow_data_merged_TN05$participant_id), i] %>%
      na.omit() %>%
      unique()
    if (length(values.tmp) == 1) {
      master_flow_data_merged_TN05[j, i] <-
        values.tmp
    } else if (length(values.tmp) > 1) {
      cat(length(values.tmp), "unique values for patient", master_flow_data_merged_TN05$participant_id[j],
          "for variable", i, "\n")
    } else if (length(values.tmp) == 0) {

      # if not found in focal data frame, look in patient_data_TN05
      if (i %in% colnames(patient_data_TN05)) {
        values.tmp <-
          patient_data_TN05[
            match(master_flow_data_merged_TN05$participant_id[j], patient_data_TN05$participant_id), i] %>%
          na.omit() %>%
          unique()
        if (length(values.tmp) == 1) {
          master_flow_data_merged_TN05[j, i] <-
            values.tmp
          next
        } else if (length(values.tmp) > 1) {
          cat(length(values.tmp), "unique values for patient", master_flow_data_merged_TN05$participant_id[j],
              "for variable", i, "\n")
          next
        } else if (length(values.tmp) == 0) {
          cat("No data for patient", master_flow_data_merged_TN05$participant_id[j],
              "for variable", i, "in patient_data_TN05\n")
        }
      } else {
        cat("Variable", i, "not found in patient_data_TN05\n")
      }

      # if not found in focal data frame or patient_data_TN05, look in cpeptide_data_TN05
      if (i %in% colnames(cpeptide_data_TN05)) {
        values.tmp <-
          cpeptide_data_TN05[
            match(master_flow_data_merged_TN05$participant_id[j], cpeptide_data_TN05$participant_id), i] %>%
          na.omit() %>%
          unique()
        if (length(values.tmp) == 1) {
          master_flow_data_merged_TN05[j, i] <-
            values.tmp
        } else if (length(values.tmp) > 1) {
          cat(length(values.tmp), "unique values for patient", master_flow_data_merged_TN05$participant_id[j],
              "for variable", i, "\n")
        } else if (length(values.tmp) == 0) {
          cat("No data for patient", master_flow_data_merged_TN05$participant_id[j],
              "for variable", i, "in cpeptide_data_TN05\n")
        }
      } else {
        cat("Variable", i, "not found in cpeptide_data_TN05\n")
      }
    }
  }
}

# fill in baseline_visit data
master_flow_data_merged_TN05$baseline_visit[is.na(master_flow_data_merged_TN05$baseline_visit)] <-
  ifelse(
    master_flow_data_merged_TN05$rnaseq_visit_number[
      is.na(master_flow_data_merged_TN05$baseline_visit)] == 2,
    "Y", "N")

# apply groups from RNAseq patients to all patients (age_split_rnaseq_breaks_all_patients)
age_split.tmp <- unique(master_flow_data_merged_TN05[,c("participant_id", "age_years")])
age_split.tmp <- age_split.tmp[age_split.tmp$age_years >= 12.1,]
age_split.tmp$age_split_rnaseq_breaks_all_patients <-
  cut(age_split.tmp$age_years, breaks=c(12.1, 19.1, 40.55), include.lowest=TRUE,
      labels=c("[12.1,19.1]", "(19.1,40.5]"))
master_flow_data_merged_TN05$age_split_rnaseq_breaks_all_patients <-
  age_split.tmp$age_split_rnaseq_breaks_all_patients[
    match(master_flow_data_merged_TN05$participant_id, age_split.tmp$participant_id)]

rm_tmp(ask=FALSE)


##### plot rate of C-peptide change vs. B cells from flow data by treatment group #####

pdf("Fig_5B.pdf", w=9, h=6)
ggplot(
  master_flow_data_merged_TN05[
    with(master_flow_data_merged_TN05,
         (baseline_visit %in% "Y")),],
  mapping=aes(
    x=cd19_pct,
    y=log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept,
    color=treatment,
    group=treatment)) +
  geom_smooth(method="lm", se=FALSE) +
  geom_point(size=3) +
  scale_color_manual(
    "Treatment",
    values=c("rituximab"=color.ritux, "placebo"=color.placebo)) +
  labs(x="% CD19+ of lymphocytes", y="Rate of C-peptide change",
       title=NULL)
dev.off()


##### fit model to rate of C-peptide change vs. CD19 percent, by treatment #####

# fit model to AUC slope by cd19_pct, by treatment, in all patients
lm.cpeptide_slope_vs_baseline_cd19_pct.by_treatment.all_patients <-
  lm(
    data=
      master_flow_data_merged_TN05[
        with(master_flow_data_merged_TN05,
             (baseline_visit %in% "Y")),],
    formula=
      log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept ~
      cd19_pct * treatment)
summary(lm.cpeptide_slope_vs_baseline_cd19_pct.by_treatment.all_patients)


##### plot rate of C-peptide change vs. B cells from flow data by treatment group, within each age split #####

# create name translators
translate_age_split2.tmp <-
  c("Younger half", "Older half") %>%
  setNames(levels(master_flow_data_merged_TN05$age_split2))

## plot cd19_pct from flow vs. rate of AUC change, split by age_split2
pdf("Fig_5D.pdf", w=9, h=6)
for (i in unique(na.omit(master_flow_data_merged_TN05$age_split2))) {
  name.tmp <- i %>%
    str_replace_all("\\[|\\]|\\(|\\)", "") %>%
    str_replace_all(",Inf", " up") %>%
    str_replace_all(",", " to ")
  print(
    ggplot(
      master_flow_data_merged_TN05[
        with(master_flow_data_merged_TN05,
             flow_visit_name=="baseline" &
               age_split2 %in% i),],
      mapping=aes(
        x=cd19_pct,
        y=log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept,
        color=treatment,
        group=treatment)) +
      geom_smooth(method="lm", se=FALSE) +
      geom_point(size=3) +
      scale_color_manual(values=c("rituximab"=color.ritux, "placebo"=color.placebo)) +
      labs(x="% CD19+ of lymphocytes", y="Rate of C-peptide change",
           title=paste0(translate_age_split2.tmp[i], " (", name.tmp, " years)")) +
      coord_cartesian(
        xlim=range(master_flow_data_merged_TN05$cd19_pct[master_flow_data_merged_TN05$rnaseq_visit_number==2], na.rm=T),
        ylim=range(
          master_flow_data_merged_TN05$log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept[
            master_flow_data_merged_TN05$rnaseq_visit_number==2], na.rm=T)))
}
dev.off()

rm_tmp(ask=FALSE)


##### fit model to rate of C-peptide change vs. CD19 percent, by treatment and age split #####

# fit model to cpeptide_slope, with cd19_pct, age_years, treatment, and interaction
lm.cpeptide_slope_vs_baseline_cd19_pct_age_split2_treatment.with_interaction <-
  lm(data=master_flow_data_merged_TN05[master_flow_data_merged_TN05$rnaseq_visit_number==2,],
     formula=
       log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept~
       cd19_pct*age_split2*treatment)
summary(lm.cpeptide_slope_vs_baseline_cd19_pct_age_split2_treatment.with_interaction)

# compute one-sided significance
contrasts.tmp <-
  matrix(c(0,0,0,0,0,1,0,0, # placebo vs. rituximab in age_split_lo
           0,0,0,0,0,1,0,1), # placebo vs. rituximab in age_split_hi
         nrow=2, byrow=TRUE)
multcomp::glht(
  lm.cpeptide_slope_vs_baseline_cd19_pct_age_split2_treatment.with_interaction,
  linfct=contrasts.tmp, alternative="greater") %>%
  summary()

contrasts.tmp <-
  matrix(c(0,0,0,0,0,0,0,1), # difference b/w (placebo vs. ritux in age_split_hi) and (placebo vs. ritux in age_split_lo)
         nrow=1, byrow=TRUE)
multcomp::glht(
  lm.cpeptide_slope_vs_baseline_cd19_pct_age_split2_treatment.with_interaction,
  linfct=contrasts.tmp, alternative="less") %>%
  summary()

rm_tmp(ask=FALSE)


##### quantify correlation of neutrophils and B cells #####

master_flow_data_merged_TN05$neutrophils_perc <-
  (sin(master_flow_data_merged_TN05$neutrophils))^2*100
master_flow_data_merged_TN05$lymphocytes_perc <-
  (sin(master_flow_data_merged_TN05$lymphocytes))^2*100

master_flow_data_merged_TN05$cd19_perc_of_cbc <-
  with(master_flow_data_merged_TN05,
       cd19_pct * lymphocytes_perc / 100)

master_flow_data_merged_TN05 %>%
  dplyr::filter(
    !is.na(neutrophils_perc),
    !is.na(cd19_perc_of_cbc),
    treatment == "placebo") %>%
  dplyr::select(neutrophils_perc, cd19_perc_of_cbc) %>%
  cor(method="pearson")

