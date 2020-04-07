##### scripts for analysis of data for Dufort et al. 2019. Cell type–specific immune phenotypes predict loss of insulin secretion in new-onset type 1 diabetes. (DOI: 10.1172/jci.insight.125556)
### this file includes scripts for analyses of complete blood count (CBC) data

##### set up environment: load packages #####

# load general packages
library(readxl)
library(colorspace)
library(RColorBrewer)
library(tidyverse)
theme_set(
  theme_bw(20) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour="black", fill=NA, size=1),
          axis.text=element_text(colour="black"),
          axis.ticks=element_line(colour="black")))
update_geom_defaults("point", list(shape=16))
library(gplots)
library(magrittr)

# load analysis_specific packages
library(edgeR)
library(limma)
library(gdata)
# library(corrMatrix)
library(lmerTest)

# load custom packages (available at github.com/benaroyaresearch)
if (!require(countSubsetNorm)) remotes::install_github("benaroyaresearch/countSubsetNorm"); library(countSubsetNorm)
if (!require(miscHelpers)) remotes::install_github("benaroyaresearch/miscHelpers"); library(miscHelpers)


##### load/save data from previous scripts #####

load(file="T1D_placebos_data_1_for_downstream_analyses.RData")


##### match merged CBC data to RNAseq and C-peptide data #####

shared_cols.tmp <- intersect(colnames(master), colnames(cbc_data))

## pull cbc_visit_number into master so that I can merge them

## use details of study to pull data into master
for (i in c("study_visit_number", "study_visit_name", "weeks",
            "cbc_visit_number", "cbc_visit_name")) {
  master[[i]] <- NA
  class(master[[i]]) <- class(study_schedules[[i]])
  for (j in 1:nrow(master)) {
    if (!is.na(master$rnaseq_visit_number[j])) {
      master[[i]][j] <-
        study_schedules[[i]][
          (study_schedules$study==master$study[j]) &
            (study_schedules$rnaseq_visit_number==master$rnaseq_visit_number[j])]
    } else if (!is.na(master$cpeptide_visit_number[j])) {
      study_schedules.tmp <-
        study_schedules[
          (study_schedules$study==master$study[j]) &
            (study_schedules$cpeptide_visit_number==master$cpeptide_visit_number[j]),]
      if (nrow(study_schedules.tmp)==1) {
        master[[i]][j] <-
          study_schedules.tmp[[i]][1]
      } else if (nrow(study_schedules.tmp) == 2) {
        master[[i]][j] <-
          study_schedules.tmp[[i]][2]
      }
    }
  }
}
rm_tmp(ask=FALSE)

## add cbc_visit_id to all rows
master$cbc_visit_id <-
  ifelse(is.na(master$cbc_visit_number), NA,
         with(master,
              paste(participant_id, cbc_visit_number, sep="_")))
any(duplicated(master$cbc_visit_id)) # no duplicates


## merge master and cbc_data
master_cbc_merged <-
  full_join(master, cbc_data)
glimpse(master_cbc_merged)

## fill in weeks from cbc_visit_number
master_cbc_merged <-
  master_cbc_merged %>%
  arrange(participant_id, cbc_visit_number)

# fill in the NAs where possible, based on cbc_visit_number
master_cbc_merged[is.na(master_cbc_merged$weeks), "weeks"] <-
  (merge(study_schedules, master_cbc_merged[is.na(master_cbc_merged$weeks),],
         all.y=TRUE, by=c("study", "cbc_visit_number"),
         sort=FALSE) %>% arrange(participant_id, cbc_visit_number))[,"weeks.x"]

## #fill in cbc_baseline_visit from cbc_visit_number

# fill in the NAs where possible, based on cbc_visit_number
master_cbc_merged[
  is.na(master_cbc_merged$cbc_baseline_visit), "cbc_baseline_visit"] <-
  with(
    master_cbc_merged[is.na(master_cbc_merged$cbc_baseline_visit),],
    ((study %in% c("ABATE", "START", "T1DAL", "TN09") & cbc_visit_number==0) |
       (study %in% c("TN02", "TN05") & cbc_visit_number==2)))

# fill in neutrophil and lymphocyte baseline values and difference (in log) from baseline
master_cbc_merged$neutrophils_abs.baseline <-
  master_cbc_merged$neutrophils_abs.logdiff <-
  master_cbc_merged$lymphocytes_abs.baseline <-
  master_cbc_merged$lymphocytes_abs.logdiff <-
  as.numeric(NA)
for (i in unique(master_cbc_merged$participant_id)) {
  data.tmp <-
    master_cbc_merged[
      with(master_cbc_merged,
           (participant_id %in% i) & cbc_baseline_visit),
      c("neutrophils_abs", "lymphocytes_abs")]
  if (nrow(data.tmp) == 1) {
    master_cbc_merged[
      master_cbc_merged$participant_id == i,
      c("neutrophils_abs.baseline", "lymphocytes_abs.baseline")] <-
      data.tmp[,c("neutrophils_abs", "lymphocytes_abs")]
  } else if (nrow(data.tmp) > 1 &
             nrow(unique(data.tmp[c("neutrophils_abs", "lymphocytes_abs")])) == 1) {
    (cat(nrow(data.tmp), "baseline rows detected for patient", i, "with identical data\n"))
    master_cbc_merged[
      master_cbc_merged$participant_id == i,
      c("neutrophils_abs.baseline", "lymphocytes_abs.baseline")] <-
      data.tmp[1,c("neutrophils_abs", "lymphocytes_abs")]
  } else (cat(nrow(data.tmp), "baseline rows detected for patient", i, "\n"))
}

master_cbc_merged[,c("neutrophils_abs.logdiff", "lymphocytes_abs.logdiff")] <-
  master_cbc_merged[,c("neutrophils_abs", "lymphocytes_abs")] -
  master_cbc_merged[,c("neutrophils_abs.baseline", "lymphocytes_abs.baseline")]

# fill in age, progressor status, and C-peptide rates (because they don't vary with time)
cols.tmp <- grep("progressor|linear|age(?!s)", colnames(master_cbc_merged), value=TRUE, perl=TRUE)
for (i in unique(master_cbc_merged$participant_id)) {
  for (j in cols.tmp) {
    data.tmp <-
      unique(na.omit(master_cbc_merged[master_cbc_merged$participant_id==i, j]))
    if (length(data.tmp) > 1) {
      if (length(unique(round(data.tmp, 3))) == 1) {
        cat("Found multiple unique values for ", j, " in Patient ", i, ", resolved by rounding.\n", sep="")
        master_cbc_merged[master_cbc_merged$participant_id==i, j] <- unique(round(data.tmp, 3))
      } else if (length(unique(round(data.tmp, 2))) == 1) {
        cat("Found multiple unique values for ", j, " in Patient ", i, ", resolved by rounding.\n", sep="")
        master_cbc_merged[master_cbc_merged$participant_id==i, j] <- unique(round(data.tmp, 2))

      } else cat("Found multiple unique values for ", j, " in Patient ", i, ", not resolved by rounding.\n", sep="")
    } else if (length(data.tmp) == 0) {
      # cat("Found no non-NA values for ", j, " in Patient ", i, ".\n", sep="")
    } else {
      master_cbc_merged[master_cbc_merged$participant_id==i, j] <- data.tmp
    }
  }
}

### scale time variables for better model fitting
## these are all scaled by the same values used in cpeptide model fitting so that they're comparable
for (var.tmp in
     grep("^(cpeptide|cbc|rnaseq)_(study|diagnosis)_day$",
          colnames(master_cbc_merged), value=TRUE)) {
  sd.tmp <- switch(str_extract(var.tmp, "day"),
                   weeks=weeks.sd,
                   day=days.sd)
  master_cbc_merged[[paste0(var.tmp, "_scaled")]] <-
    scale(master_cbc_merged[[var.tmp]], center=FALSE,
          scale=sd.tmp)[,1]
}

rm_tmp(ask=FALSE)


##### fit model of percent neutrophils to fast/slow progressor status #####

## neutrophil percent by cbc_study_day and progressor status
lmer.perc_neutrophils.cbc_study_day.progressor_split_cpep_rate_fast25_slow75 <-
  lmer(
    neutrophils ~
      cbc_study_day + progressor_split_cpep_rate_fast25_slow75 + (1|participant_id),
    data=
      master_cbc_merged[!duplicated(master_cbc_merged$cbc_visit_id),])
summary(lmer.perc_neutrophils.cbc_study_day.progressor_split_cpep_rate_fast25_slow75)

## neutrophil percent by cbc_study_day and progressor status, with interaction
lmer.perc_neutrophils.cbc_study_day.progressor_split_cpep_rate_fast25_slow75.inc_interaction <-
  lmer(
    neutrophils ~
      cbc_study_day * progressor_split_cpep_rate_fast25_slow75 + (1|participant_id),
    data=master_cbc_merged[!duplicated(master_cbc_merged$cbc_visit_id),])
summary(lmer.perc_neutrophils.cbc_study_day.progressor_split_cpep_rate_fast25_slow75.inc_interaction)


##### plot % neutrophils over time by fast/slow progressor status #####

pdf("Fig_3A.pdf", w=10,h=6)
ggplot(
  data=
    master_cbc_merged[
      !duplicated(master_cbc_merged$cbc_visit_id) &
        !is.na(master_cbc_merged$progressor_split_cpep_rate_fast25_slow75) &
        !is.na(master_cbc_merged$neutrophils),] %>%
    cbind(
      neutrophils_predicted =
        predict(lmer.perc_neutrophils.cbc_study_day.progressor_split_cpep_rate_fast25_slow75)),
  mapping=aes(
    x=cbc_study_day, y=neutrophils,
    group=progressor_split_cpep_rate_fast25_slow75,
    colour=progressor_split_cpep_rate_fast25_slow75)) +
  geom_point(alpha=0.5, size=3, shape=16) +
  geom_smooth(
    method="lm",
    mapping=aes(y=neutrophils_predicted),
    size=2.0, se=FALSE, color="black") +
  geom_smooth(
    method="lm",
    mapping=aes(y=neutrophils_predicted),
    size=1.5, se=FALSE) +
  scale_color_manual(name="Progressor", values=c("red", "blue")) +
  scale_x_continuous(breaks=seq(0,1500, by=500)) +
  labs(y="Neutrophil %\n(arcsine-transformed)", x="Days in study")
dev.off()


##### fit model of percent lymphocytes to fast/slow progressor status #####

## lymphocyte percent by cbc_study_day and progressor status
lmer.perc_lymphocytes.cbc_study_day.progressor_split_cpep_rate_fast25_slow75 <-
  lmer(
    lymphocytes ~
      cbc_study_day + progressor_split_cpep_rate_fast25_slow75 + (1|participant_id),
    data=master_cbc_merged[!duplicated(master_cbc_merged$cbc_visit_id),])
summary(lmer.perc_lymphocytes.cbc_study_day.progressor_split_cpep_rate_fast25_slow75)

## neutrophil percent by cbc_study_day and progressor status, with interaction
lmer.perc_lymphocytes.cbc_study_day.progressor_split_cpep_rate_fast25_slow75.inc_interaction <-
  lmer(
    lymphocytes ~
      cbc_study_day * progressor_split_cpep_rate_fast25_slow75 + (1|participant_id),
    data=master_cbc_merged[!duplicated(master_cbc_merged$cbc_visit_id),])
summary(lmer.perc_lymphocytes.cbc_study_day.progressor_split_cpep_rate_fast25_slow75.inc_interaction)


##### plot % lymphocytes over time by fast/slow progressor status #####

pdf("Fig_S9A.pdf", w=10, h=6)
ggplot(
  data=
    master_cbc_merged[
      !duplicated(master_cbc_merged$cbc_visit_id) &
        !is.na(master_cbc_merged$progressor_split_cpep_rate_fast25_slow75) &
        !is.na(master_cbc_merged$lymphocytes),] %>%
    cbind(
      lymphocytes_predicted =
        predict(lmer.perc_lymphocytes.cbc_study_day.progressor_split_cpep_rate_fast25_slow75)),
  mapping=aes(
    x=cbc_study_day, y=lymphocytes,
    group=progressor_split_cpep_rate_fast25_slow75,
    colour=progressor_split_cpep_rate_fast25_slow75)) +
  geom_point(alpha=0.5, size=2, shape=16) +
  geom_smooth(
    method="lm",
    mapping=aes(y=lymphocytes_predicted),
    size=2.0, se=FALSE, color="black") +
  geom_smooth(
    method="lm",
    mapping=aes(y=lymphocytes_predicted),
    size=1.5, se=FALSE) +
  scale_color_manual(name="Progressor", values=c("red", "blue")) +
  scale_x_continuous(breaks=seq(0,1500, by=500)) +
  labs(y="Lymphocyte %\n(arcsine-transformed)", x="Days in study")
dev.off()


##### fit model of lymphocyte-to-neutrophil ratio to fast/slow progressor status #####

## lymphocytes_neutrophils_ratio by cbc_study_day and progressor status
lmer.lymphocytes_neutrophils_ratio.cbc_study_day.progressor_split_cpep_rate_fast25_slow75 <-
  lmer(
    lymphocytes_neutrophils_ratio ~
      cbc_study_day + progressor_split_cpep_rate_fast25_slow75 + (1|participant_id),
    data=master_cbc_merged[!duplicated(master_cbc_merged$cbc_visit_id),])
summary(lmer.lymphocytes_neutrophils_ratio.cbc_study_day.progressor_split_cpep_rate_fast25_slow75)

## lymphocytes_neutrophils_ratio by cbc_study_day and progressor status
lmer.lymphocytes_neutrophils_ratio.cbc_study_day.progressor_split_cpep_rate_fast25_slow75.inc_interaction <-
  lmer(
    lymphocytes_neutrophils_ratio ~
      cbc_study_day * progressor_split_cpep_rate_fast25_slow75 + (1|participant_id),
    data=master_cbc_merged[!duplicated(master_cbc_merged$cbc_visit_id),])
summary(lmer.lymphocytes_neutrophils_ratio.cbc_study_day.progressor_split_cpep_rate_fast25_slow75.inc_interaction)


##### plot lymphocyte-to-neutrophil ratio over time by fast/slow progressor status #####

data.tmp <-
  master_cbc_merged[
    !duplicated(master_cbc_merged$cbc_visit_id) &
      !is.na(master_cbc_merged$progressor_split_cpep_rate_fast25_slow75) &
      !is.na(master_cbc_merged$lymphocytes_neutrophils_ratio),] %>%
  cbind(
    lymphocytes_neutrophils_ratio_predicted =
      predict(lmer.lymphocytes_neutrophils_ratio.cbc_study_day.progressor_split_cpep_rate_fast25_slow75))

pdf("Fig_S9B.pdf", w=10, h=6)
ggplot(
  data=
    master_cbc_merged[
      !duplicated(master_cbc_merged$cbc_visit_id) &
        !is.na(master_cbc_merged$progressor_split_cpep_rate_fast25_slow75) &
        !is.na(master_cbc_merged$lymphocytes_neutrophils_ratio),] %>%
    cbind(
      lymphocytes_neutrophils_ratio_predicted =
        predict(lmer.lymphocytes_neutrophils_ratio.cbc_study_day.progressor_split_cpep_rate_fast25_slow75)),
  mapping=aes(
    x=cbc_study_day, y=lymphocytes_neutrophils_ratio,
    group=progressor_split_cpep_rate_fast25_slow75,
    colour=progressor_split_cpep_rate_fast25_slow75)) +
  geom_point(alpha=0.5, size=2, shape=16) +
  geom_smooth(
    method="lm",
    mapping=aes(y=lymphocytes_neutrophils_ratio_predicted),
    size=2.0, se=FALSE, color="black") +
  geom_smooth(
    method="lm",
    mapping=aes(y=lymphocytes_neutrophils_ratio_predicted),
    size=1.5, se=FALSE) +
  scale_color_manual(name="Progressor", values=c("red", "blue")) +
  scale_x_continuous(breaks=seq(0,1500, by=500)) +
  labs(y="Lymphocyte:Neutrophil Ratio", x="Days in study")
dev.off()


##### remove excess time columns from master_cbc_merged #####

master_cbc_merged <-
  master_cbc_merged[,!str_detect(colnames(master_cbc_merged), "scaled_")]


##### combine cell subset models into a single object for easier storage #####

cbc_models.all <- list()
for (i in ls_grep("^lm")) cbc_models.all[[i]] <- get(i)

rm_tmp(ask=FALSE)


##### filter and transform counts object and master_cbc_merged object #####

counts.filtered_cbcs <-
  calc_norm_counts(
    counts=counts.final, design=master_cbc_merged, libID_col="libid",
    min_cpm=0, min_libs_perc = 0,
    normalize=FALSE,
    return_DGEcounts=FALSE)
master_cbc_merged.filtered_cbcs <-
  master_cbc_merged[match(colnames(counts.filtered_cbcs), master_cbc_merged$libid),]


##### save objects for downstream use #####

save(file="T1D_placebos_data_2_for_downstream_analyses.RData",
     list=c(
       "master_cbc_merged", "master_cbc_merged.filtered_cbcs", "counts.filtered_cbcs",
       "counts.final", "cbc_models.all", "vwts.all"))
