## scripts to analyze data from the TN-05 rituximab trial, as part of the T1D placebos projet

##### set up environment: load packages #####

## load general packages
library(xlsx)
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
library(ggbeeswarm)

## load analysis-specific packages
library(limma)
library(edgeR)
library(lmerTest)
library(annotables)
library(NOISeq)
library(datamart)
library(ordinal)

# load useful functions
library(RNAseQC)
library(countSubsetNorm)
library(miscHelpers)
library(geneSetTools)
library(limmaTools)


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
patient_data.TN05 <- read.csv("TN05_patient_data.csv")

## CBC data
cbc.TN05 <- read.csv("TN05_cbc_data_processed_merged.csv")

## C-peptide data
cpeptide_data.TN05 <- read.csv("TN05_cpeptide_data_processed_merged.csv")


##### calculate C-peptide derivative values #####

## determine baseline values, and calculate log auc
cols.tmp <-
  c("baseline_auc2hr", "baseline_log_auc2hr")
cpeptide_data.TN05[, cols.tmp] <- as.numeric(NA)

for (i in unique(cpeptide_data.TN05$participant_id)) {
  rows.tmp <- which(cpeptide_data.TN05$participant_id==i)
  data.tmp <- cpeptide_data.TN05[rows.tmp,]
  baseline_row.tmp <- which(data.tmp$cpeptide_visit_name %in% c("baseline", "0"))
  
  if (length(baseline_row.tmp)==0)
    baseline_row.tmp <- which(data.tmp$cpeptide_visit_name %in% c("screening", "-1"))
  
  if (length(baseline_row.tmp)==0) next
  
  cpeptide_data.TN05[rows.tmp, cols.tmp] <-
    data.tmp[baseline_row.tmp, str_replace(cols.tmp, "baseline_", "")]
}

rm_tmp(ask=FALSE)


##### determine visits with C-peptide previously at detection threshold, to exclude them #####

low_lim_detect <- 3 * ng_ml_to_nmol_l_min

cpeptide_data.TN05$prev_at_cpep_detect_thresh <- as.logical(NA)
cpeptide_data.TN05 <-
  cpeptide_data.TN05 %>%
  arrange(participant_id, cpeptide_study_day)

for (i in unique(cpeptide_data.TN05$participant_id)) { # iterate over patient
  rows.tmp <-
    which(
      with(cpeptide_data.TN05,
           participant_id==i & !is.na(auc2hr))) # determine all rows for current patient
  cpeptide_data.TN05$prev_at_cpep_detect_thresh[rows.tmp[1]] <- FALSE # set first visit to FALSE
  if (length(rows.tmp) > 1)
    for (j in (2:length(rows.tmp))) {
      cpeptide_data.TN05$prev_at_cpep_detect_thresh[rows.tmp[j]] <-
        (cpeptide_data.TN05$auc2hr[rows.tmp[j-1]] < low_lim_detect * 1.01) |  # previous visit is at detection limit
        cpeptide_data.TN05$prev_at_cpep_detect_thresh[rows.tmp[j-1]]  # previous visit had prev_at_cpep_detect_thresh
    }
}

rm_tmp(ask=FALSE)


##### fit curves to c-peptide AUC by cpeptide_study_day #####

# scale cpeptide_study_day so that it works better in model fitting
cpeptide_study_day_scaled.tmp <-
  scale(cpeptide_data.TN05$cpeptide_study_day, center=FALSE,
        scale=sd(cpeptide_data.TN05$cpeptide_study_day, na.rm=TRUE))
cpeptide_data.TN05$cpeptide_study_day_scaled <-
  as.vector(cpeptide_study_day_scaled.tmp)
days.sd <- attr(cpeptide_study_day_scaled.tmp, "scaled:scale")

## fit a mixed-effects model with a term to allow different fixed-effect slopes by treatment
lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random <-
  lmer(
    log_auc2hr ~
      (1|participant_id) +
      cpeptide_study_day_scaled +
      cpeptide_study_day_scaled:treatment +
      (cpeptide_study_day_scaled-1|participant_id),
    data=cpeptide_data.TN05[!cpeptide_data.TN05$prev_at_cpep_detect_thresh,])
summary(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)


##### extract slopes from model fit to log_auc2hr by cpeptide_study_day #####

cpeptide_data.TN05$log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept <-
  as.numeric(NA)

for (i in 1:nrow(cpeptide_data.TN05)) {
  # extract slope for linear random with intercept, with separate fixed effects for placebo and rituximab
  if (sum(!is.na(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id[
    match(
      cpeptide_data.TN05$participant_id[i],
      rownames(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id)),])) == 3)
    cpeptide_data.TN05$log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept[i] <-
      (coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id$cpeptide_study_day_scaled[
        match(cpeptide_data.TN05$participant_id[i],
              rownames(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id))] +
         coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id$'cpeptide_study_day_scaled:treatmentrituximab'[
           match(cpeptide_data.TN05$participant_id[i],
                 rownames(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id))] *
         (cpeptide_data.TN05$treatment[i]=="rituximab")) /
      days.sd * 365.25
}

rm_tmp(ask=FALSE)


##### plot and model rate of C-peptide change vs. age split and treatment, with interaction #####

## plot C-peptide rate
pdf("Fig_S12.pdf", w=11, h=6)
ggplot(
  cpeptide_data.TN05[!duplicated(cpeptide_data.TN05$participant_id),],
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
  lm(data=cpeptide_data.TN05[!duplicated(cpeptide_data.TN05$participant_id),],
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

# read in filtered, normalized counts
counts.TN05 <-
  read.table(
    "TN05_counts_normalized.txt",
    sep="\t", header=T)

## read in sample annotation
sample_annotation.TN05 <- read.csv("TN05_sample_annotation.csv")

# merge patient_data.TN05 into sample_annotation.TN05 data
sample_annotation.TN05 <-
  merge(sample_annotation.TN05,
        patient_data.TN05[
          c("participant_id",
            setdiff(colnames(patient_data.TN05), colnames(sample_annotation.TN05)))],
        by="participant_id", all.x=T)

# generate visit_id 
sample_annotation.TN05$rnaseq_visit_id <-
  with(sample_annotation.TN05,
       paste(participant_id, rnaseq_visit_number, sep="_"))
sample_annotation.TN05$cpeptide_visit_id <-
  with(sample_annotation.TN05,
       paste(participant_id, cpeptide_visit_number, sep="_"))

# merge cpeptide_data.TN05 into sample_annotation.TN05
sample_annotation.TN05 <-
  merge(sample_annotation.TN05,
        cpeptide_data.TN05[
          c("cpeptide_visit_id",
            setdiff(colnames(cpeptide_data.TN05), colnames(sample_annotation.TN05)))],
        by="cpeptide_visit_id", all.x=T)

rm_tmp(ask=FALSE)


##### filter sample_annotation.TN05 and counts.TN05 to include only shared libraries, and put them in the same order #####

sample_annotation.TN05.merged <-
  sample_annotation.TN05[sample_annotation.TN05$libid %in% colnames(counts.TN05),] %>%
  arrange(participant_id, rnaseq_visit_number)
counts.TN05.merged <- counts.TN05[,match(sample_annotation.TN05.merged$libid, colnames(counts.TN05))]


##### some basic checks for problematic data #####

# calculate ratio of reads mapping to X and to Y chromosome, for each library
logXY.tmp <- logXYratio(counts.TN05.merged, lib_cols=1:ncol(counts.TN05.merged), gene_ID="symbol")
hist(logXY.tmp, breaks=15, main="log-ratio of X reads to Y reads")
cbind(sample_annotation.TN05.merged$sex,
      logXY.tmp[match(sample_annotation.TN05.merged$libid, names(logXY.tmp))])[
        order(logXY.tmp[match(sample_annotation.TN05.merged$libid, names(logXY.tmp))]),]
sort(logXY.tmp)
logXY_threshold <- 4

sample_annotation.TN05.merged$sex_by_rna <-
  ifelse(logXY.tmp[sample_annotation.TN05.merged$libid] > logXY_threshold, "F", "M")

# check that inferred sex is same for all libraries for each patient
table(sample_annotation.TN05.merged[,c("participant_id", "sex_by_rna")])
which(rowSums(table(sample_annotation.TN05.merged[,c("participant_id", "sex_by_rna")])==0)!=1)
# looks good!

# check that inferred sex matches sex from database
all(sample_annotation.TN05.merged$sex == sample_annotation.TN05.merged$sex_by_rna)
# looks good!

# drop problematic libraries, and create quality-controlled sample_annotation.TN05 and counts.TN05 objects
# except no problematic libraries to drop here
sample_annotation.TN05.merged.qc <- sample_annotation.TN05.merged
counts.TN05.merged.qc <- counts.TN05.merged[,match(sample_annotation.TN05.merged.qc$libid, colnames(counts.TN05.merged))]

rm_tmp(ask=FALSE)


##### generate simply-named combined objects with all libraries passing QC and outlier cuts #####

counts.final.TN05 <- counts.TN05.merged.qc
master.TN05 <- sample_annotation.TN05.merged.qc

# split subjects by age above/below median
age_split.tmp <- unique(master.TN05[,c("participant_id", "age_years")])
age_split.tmp$age_split_rnaseq_only <- cut_number(age_split.tmp$age_years, n=2)
master.TN05$age_split_rnaseq_only <-
  age_split.tmp$age_split_rnaseq_only[
    match(master.TN05$participant_id, age_split.tmp$participant_id)]
master.TN05[,c("age_split_rnaseq_only", "age_years")] # looks good!

### merge CBC data into master.TN05
master.TN05 <-
  base::merge(
    master.TN05,
    cbc.TN05[,setdiff(colnames(cbc.TN05), colnames(master.TN05))],
    by.x="rnaseq_visit_id",
    by.y="cbc_visit_id",
    all.x=TRUE) %>%
  arrange(rnaseq_visit_id)


##### calculate median expression of CD19.mod #####

vwts.all.TN05 <-
  calc_norm_counts(
    counts=counts.final.TN05, design=master.TN05, libID_col="libid",
    min_cpm=1, min_libs_perc=0.15,
    normalize=FALSE,
    return_DGEcounts=TRUE) %>%
  voom()

## store median gene set expression in master.TN05
master.TN05$median_CD19.mod <-
  gene_set_median_count(
    gene_sets.Linsley[["CD19.mod"]], vwts.all.TN05,
    remove_low_count_genes=TRUE)[
      match(master.TN05$libid, colnames(vwts.all.TN05$E))]

## determine baseline median scores for CD19.mod
master.TN05[, "baseline_median_CD19.mod"] <- as.numeric(NA)
for (i in unique(master.TN05$participant_id)) {
  master.tmp <-
    master.TN05[
      with(master.TN05,
           (participant_id==i) & (rnaseq_visit_number==2)),]
  if (nrow(master.tmp)==1) {
    master.TN05[
      master.TN05$participant_id==i,
      "baseline_median_CD19.mod"] <-
      master.tmp[, "median_CD19.mod"]
  } else cat(nrow(master.tmp), "baseline visits for", i, "\n")
}

##### plot and model rate of C-peptide change vs. median gene set expression, by treatment group #####

## plot it
pdf("Fig_5A.pdf", w=9, h=6)
print(
  ggplot(
    master.TN05[
      with(master.TN05,
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
      xlim=range(master.TN05$median_CD19.mod[master.TN05$rnaseq_visit_number==2], na.rm=T),
      ylim=range(
        master.TN05$log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept[
          master.TN05$rnaseq_visit_number==2], na.rm=T)))
dev.off()

## fit model to AUC slope by median_CD19.mod, by treatment
lm.cpeptide_slope_vs_baseline_median_CD19.mod.by_treatment <-
  lm(data=master.TN05[master.TN05$rnaseq_visit_number==2,],
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
for (i in unique(master.TN05$age_split_rnaseq_only)) {
  name.tmp <- i %>%
    str_replace_all("\\[|\\]|\\(|\\)", "") %>%
    str_replace_all(",Inf", " up") %>%
    str_replace_all(",", " to ")
  print(
    ggplot(
      master.TN05[
        with(master.TN05,
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
        xlim=range(master.TN05$median_CD19.mod[master.TN05$rnaseq_visit_number==2], na.rm=T),
        ylim=range(
          master.TN05$log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept[
            master.TN05$rnaseq_visit_number==2], na.rm=T)))
}
dev.off()

## fit model to cpeptide_slope, with median_CD19.mod, age_split_rnaseq_only, treatment, and interaction
lm.cpeptide_slope_vs_baseline_median_CD19.mod_age_split_rnaseq_only_treatment.with_interaction <-
  lm(data=master.TN05[master.TN05$rnaseq_visit_number==2,],
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

flow_data.TN05 <- read.csv("TN05_flow_data_processed_merged.csv")

## merge patient_data.TN05 into flow_data.TN05
flow_data.TN05 <-
  merge(flow_data.TN05, patient_data.TN05,
        all.x=TRUE)

## populate additional data fields in master.flow_data.merged
flow_data.TN05$flow_study_day <-
  with(flow_data.TN05,
       flow_collection_date - day_0)

flow_data.TN05$baseline_visit <-
  ifelse(flow_data.TN05$flow_visit_number==2, "Y", "N")

## merge cpeptide_data.TN05 into flow_data.TN05
flow_cpeptide_data.merged <-
  merge(flow_data.TN05,
        cpeptide_data.TN05,
        all.x=TRUE)


##### merge flow_data with master RNAseq annotation object #####

intersect(colnames(flow_data.TN05), colnames(master.TN05))
master.flow_data.merged.TN05 <-
  merge(master.TN05,
        flow_data.TN05[
          ,c("study", "random_id", "participant_id",
             setdiff(colnames(flow_data.TN05), colnames(master.TN05)))],
        by.x=c("study", "random_id", "participant_id", "rnaseq_visit_number"),
        by.y=c("study", "random_id", "participant_id", "study_visit_number"),
        all=TRUE)

# check for all-NA rows
ncol(master.flow_data.merged.TN05)
rowSums(!is.na(master.flow_data.merged.TN05)) %>%
  min()
# none; one with 21 non-NA columns

# propagate patient-specific values (e.g. progressor status and C-peptide slope) to all datapoints for each patient
cols.to_propagate.tmp <-
  c("participant_id", "group", "treatment", "age_years", "gender",
    "race", "ethnicity", "sex", "age_split2", "age_split_rnaseq_only",
    "baseline_auc2hr", "baseline_log_auc2hr",
    "log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept")

colSums(is.na(master.flow_data.merged.TN05))
colSums(is.na(master.flow_data.merged.TN05[,cols.to_propagate.tmp])) %>%
  sum()

for (i in cols.to_propagate.tmp) {
  rows.tmp <- which(is.na(master.flow_data.merged.TN05[,i]))
  cat("Rows to fill for variable", i, ": ", rows.tmp, "\n")
  for (j in rows.tmp) {
    # first try to get data from within data frame
    values.tmp <-
      master.flow_data.merged.TN05[
        match(master.flow_data.merged.TN05$random_id[j], master.flow_data.merged.TN05$random_id), i] %>%
      na.omit() %>%
      unique()
    if (length(values.tmp) == 1) {
      master.flow_data.merged.TN05[j, i] <-
        values.tmp
    } else if (length(values.tmp) > 1) {
      cat(length(values.tmp), "unique values for patient", master.flow_data.merged.TN05$random_id[j],
          "for variable", i, "\n")
    } else if (length(values.tmp) == 0) {
      
      # if not found in focal data frame, look in patient_data.TN05
      if (i %in% colnames(patient_data.TN05)) {
        values.tmp <- 
          patient_data.TN05[
            match(master.flow_data.merged.TN05$random_id[j], patient_data.TN05$random_id), i] %>%
          na.omit() %>%
          unique()
        if (length(values.tmp) == 1) {
          master.flow_data.merged.TN05[j, i] <-
            values.tmp
          next
        } else if (length(values.tmp) > 1) {
          cat(length(values.tmp), "unique values for patient", master.flow_data.merged.TN05$random_id[j],
              "for variable", i, "\n")
          next
        } else if (length(values.tmp) == 0) {
          cat("No data for patient", master.flow_data.merged.TN05$random_id[j],
              "for variable", i, "in patient_data.TN05\n")
        }
      } else {
        cat("Variable", i, "not found in patient_data.TN05\n")
      }
      
      # if not found in focal data frame or patient_data.TN05, look in cpeptide_data.TN05
      if (i %in% colnames(cpeptide_data.TN05)) {
        values.tmp <- 
          cpeptide_data.TN05[
            match(master.flow_data.merged.TN05$random_id[j], cpeptide_data.TN05$random_id), i] %>%
          na.omit() %>%
          unique()
        if (length(values.tmp) == 1) {
          master.flow_data.merged.TN05[j, i] <-
            values.tmp
        } else if (length(values.tmp) > 1) {
          cat(length(values.tmp), "unique values for patient", master.flow_data.merged.TN05$random_id[j],
              "for variable", i, "\n")
        } else if (length(values.tmp) == 0) {
          cat("No data for patient", master.flow_data.merged.TN05$random_id[j],
              "for variable", i, "in cpeptide_data.TN05\n")
        }
      } else {
        cat("Variable", i, "not found in cpeptide_data.TN05\n")
      }
    }
  }
}

# fill in baseline_visit data
master.flow_data.merged.TN05$baseline_visit[is.na(master.flow_data.merged.TN05$baseline_visit)] <-
  ifelse(
    master.flow_data.merged.TN05$rnaseq_visit_number[
      is.na(master.flow_data.merged.TN05$baseline_visit)] == 2,
    "Y", "N")

# apply groups from RNAseq patients to all patients (age_split_rnaseq_breaks_all_patients)
age_split.tmp <- unique(master.flow_data.merged.TN05[,c("participant_id", "age_years")])
age_split.tmp <- age_split.tmp[age_split.tmp$age_years >= 12.1,]
age_split.tmp$age_split_rnaseq_breaks_all_patients <-
  cut(age_split.tmp$age_years, breaks=c(12.1, 19.1, 40.55), include.lowest=TRUE,
      labels=c("[12.1,19.1]", "(19.1,40.5]"))
master.flow_data.merged.TN05$age_split_rnaseq_breaks_all_patients <-
  age_split.tmp$age_split_rnaseq_breaks_all_patients[
    match(master.flow_data.merged.TN05$participant_id, age_split.tmp$participant_id)]

rm_tmp(ask=FALSE)


##### plot rate of C-peptide change vs. B cells from flow data by treatment group #####

pdf("Fig_5B.pdf", w=9, h=6)
ggplot(
  master.flow_data.merged.TN05[
    with(master.flow_data.merged.TN05,
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
      master.flow_data.merged.TN05[
        with(master.flow_data.merged.TN05,
             (baseline_visit %in% "Y")),],
    formula=
      log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept ~
      cd19_pct * treatment)
summary(lm.cpeptide_slope_vs_baseline_cd19_pct.by_treatment.all_patients)


##### plot rate of C-peptide change vs. B cells from flow data by treatment group, within each age split #####

# create name translators
translate_age_split2.tmp <-
  c("Younger half", "Older half") %>%
  setNames(levels(master.flow_data.merged.TN05$age_split2))

## plot cd19_pct from flow vs. rate of AUC change, split by age_split2
pdf("Fig_5D.pdf", w=9, h=6)
for (i in unique(na.omit(master.flow_data.merged.TN05$age_split2))) {
  name.tmp <- i %>%
    str_replace_all("\\[|\\]|\\(|\\)", "") %>%
    str_replace_all(",Inf", " up") %>%
    str_replace_all(",", " to ")
  print(
    ggplot(
      master.flow_data.merged.TN05[
        with(master.flow_data.merged.TN05,
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
           title=paste0(translate_age_split2[i], " (", name.tmp, " years)")) +
      coord_cartesian(
        xlim=range(master.flow_data.merged.TN05$cd19_pct[master.flow_data.merged.TN05$rnaseq_visit_number==2], na.rm=T),
        ylim=range(
          master.flow_data.merged.TN05$log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept[
            master.flow_data.merged.TN05$rnaseq_visit_number==2], na.rm=T)))
}
dev.off()

rm_tmp(ask=FALSE)


##### fit model to rate of C-peptide change vs. CD19 percent, by treatment and age split #####

# fit model to cpeptide_slope, with cd19_pct, age_years, treatment, and interaction
lm.cpeptide_slope_vs_baseline_cd19_pct_age_split2_treatment.with_interaction <-
  lm(data=master.flow_data.merged.TN05[master.flow_data.merged.TN05$rnaseq_visit_number==2,],
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

master.flow_data.merged.TN05$neutrophils_perc <-
  (sin(master.flow_data.merged.TN05$neutrophils))^2*100
master.flow_data.merged.TN05$lymphocytes_perc <-
  (sin(master.flow_data.merged.TN05$lymphocytes))^2*100

master.flow_data.merged.TN05$cd19_perc_of_cbc <-
  with(master.flow_data.merged.TN05,
       cd19_pct * lymphocytes_perc / 100)

master.flow_data.merged.TN05 %>%
  dplyr::filter(
    !is.na(neutrophils_perc),
    !is.na(cd19_perc_of_cbc),
    treatment == "placebo") %>%
  dplyr::select(neutrophils_perc, cd19_perc_of_cbc) %>%
  cor(method="pearson")
