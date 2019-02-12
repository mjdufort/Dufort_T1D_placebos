## scripts to analyze cell blood count (CBC) data for the T1D placebos project

##### set up environment: load packages #####

# load general packages
library(xlsx)
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
library(corrMatrix)
library(lmerTest)

# load my relevant functions
library(countSubsetNorm)
library(miscHelpers)


##### load/save data from previous scripts #####

# load data from previous CBC data processing
load(file="T1D_placebos_data_1_for_downstream_analyses.RData")


##### match merged CBC data to RNAseq and C-peptide data #####

shared_cols.tmp <- intersect(colnames(master), colnames(cbc.merged))
master.cbc.merged <-
  base::merge(
    master, cbc.merged, all=TRUE,
    by=c("cbc_visit_id", "cbc_visit_date"))
glimpse(master.cbc.merged)

# fill in data where missing in one shared column
for (i in grep("\\.x", colnames(master.cbc.merged), value=TRUE)) {
  master.cbc.merged[is.na(master.cbc.merged[[i]]), i] <-
    master.cbc.merged[is.na(master.cbc.merged[[i]]), str_replace(i, "\\.x", ".y")]
}

# drop and rename shared columns
master.cbc.merged <-
  master.cbc.merged[, !str_detect(colnames(master.cbc.merged), "\\.y")]
colnames(master.cbc.merged) <-
  colnames(master.cbc.merged) %>%
  str_replace("\\.x", "")

## fill in weeks from cbc_visit_number
master.cbc.merged <-
  master.cbc.merged %>%
  arrange(participant_id, cbc_visit_number)

# fill in the NAs where possible, based on cbc_visit_number
master.cbc.merged[is.na(master.cbc.merged$weeks), "weeks"] <-
  (merge(study_schedules, master.cbc.merged[is.na(master.cbc.merged$weeks),],
         all.y=TRUE, by=c("study", "cbc_visit_number"),
         sort=FALSE) %>% arrange(participant_id, cbc_visit_number))[,"weeks.x"]

## #fill in cbc_baseline_visit from cbc_visit_number

# fill in the NAs where possible, based on cbc_visit_number
master.cbc.merged[
  is.na(master.cbc.merged$cbc_baseline_visit), "cbc_baseline_visit"] <-
  with(
    master.cbc.merged[is.na(master.cbc.merged$cbc_baseline_visit),],
    ((study %in% c("ABATE", "START", "T1DAL", "TN09") & cbc_visit_number==0) |
       (study %in% c("TN02", "TN05") & cbc_visit_number==2)))

# fill in neutrophil and lymphocyte baseline values and difference (in log) from baseline
master.cbc.merged$neutrophils_abs.baseline <- 
  master.cbc.merged$neutrophils_abs.logdiff <- 
  master.cbc.merged$lymphocytes_abs.baseline <-
  master.cbc.merged$lymphocytes_abs.logdiff <-
  as.numeric(NA)
for (i in unique(master.cbc.merged$participant_id)) {
  data.tmp <-
    master.cbc.merged[
      with(master.cbc.merged,
           (participant_id %in% i) & cbc_baseline_visit),
      c("neutrophils_abs", "lymphocytes_abs")]
  if (nrow(data.tmp) == 1) {
    master.cbc.merged[
      master.cbc.merged$participant_id == i,
      c("neutrophils_abs.baseline", "lymphocytes_abs.baseline")] <-
      data.tmp[,c("neutrophils_abs", "lymphocytes_abs")]
  } else if (nrow(data.tmp) > 1 &
             nrow(unique(data.tmp[c("neutrophils_abs", "lymphocytes_abs")])) == 1) {
    (cat(nrow(data.tmp), "baseline rows detected for patient", i, "with identical data\n"))
    master.cbc.merged[
      master.cbc.merged$participant_id == i,
      c("neutrophils_abs.baseline", "lymphocytes_abs.baseline")] <-
      data.tmp[1,c("neutrophils_abs", "lymphocytes_abs")]
  } else (cat(nrow(data.tmp), "baseline rows detected for patient", i, "\n"))
}

master.cbc.merged[,c("neutrophils_abs.logdiff", "lymphocytes_abs.logdiff")] <-
  master.cbc.merged[,c("neutrophils_abs", "lymphocytes_abs")] -
  master.cbc.merged[,c("neutrophils_abs.baseline", "lymphocytes_abs.baseline")]

# fill in age, progressor status, and C-peptide rates (because they don't vary with time)
cols.tmp <- grep("progressor|linear|age(?!s)", colnames(master.cbc.merged), value=TRUE, perl=TRUE)
for (i in unique(master.cbc.merged$participant_id)) {
  for (j in cols.tmp) {
    data.tmp <-
      unique(na.omit(master.cbc.merged[master.cbc.merged$participant_id==i, j]))
    if (length(data.tmp) > 1) {
      if (length(unique(round(data.tmp, 3))) == 1) {
        cat("Found multiple unique values for ", j, " in Patient ", i, ", resolved by rounding.\n", sep="")
        master.cbc.merged[master.cbc.merged$participant_id==i, j] <- unique(round(data.tmp, 3))
      } else if (length(unique(round(data.tmp, 2))) == 1) {
        cat("Found multiple unique values for ", j, " in Patient ", i, ", resolved by rounding.\n", sep="")
        master.cbc.merged[master.cbc.merged$participant_id==i, j] <- unique(round(data.tmp, 2))
        
      } else cat("Found multiple unique values for ", j, " in Patient ", i, ", not resolved by rounding.\n", sep="")
    } else if (length(data.tmp) == 0) {
      # cat("Found no non-NA values for ", j, " in Patient ", i, ".\n", sep="")
    } else {
      master.cbc.merged[master.cbc.merged$participant_id==i, j] <- data.tmp
    }
  }
}

### scale time variables for better model fitting
## these are all scaled by the comparable scales used in cpeptide model fitting so that they're comparable
for (var.tmp in
     grep("^(cpeptide|cbc|rnaseq)_(study|diagnosis)_day$",
          colnames(master.cbc.merged), value=TRUE)) {
  sd.tmp <- switch(str_extract(var.tmp, "day"),
                   weeks=weeks.sd,
                   day=days.sd)
  master.cbc.merged[[paste0(var.tmp, "_scaled")]] <-
    scale(master.cbc.merged[[var.tmp]], center=FALSE,
          scale=sd.tmp)[,1]
}

rm_tmp(ask=FALSE)


##### check overlap in data types (RNAseq, CBCs, C-peptide AUC) #####

# counts of RNAseq data by study and timepoint
table(master.final[!is.na(master.final$libid),c("study","weeks")])

# counts of C-peptide AUC + RNAseq data by study and timepoint
table(master[!is.na(master$libid),c("study","weeks")])


##### fit model of percent neutrophils to fast/slow progressor status #####

## neutrophil percent by cbc_study_day and progressor status
lmer.perc_neutrophils.cbc_study_day.progressor_split_cpep_rate_fast25_slow75 <-
  lmer(
    neutrophils ~
      cbc_study_day + progressor_split_cpep_rate_fast25_slow75 + (1|participant_id),
    data=master.cbc.merged[!duplicated(master.cbc.merged$cbc_visit_id_date),])
summary(lmer.perc_neutrophils.cbc_study_day.progressor_split_cpep_rate_fast25_slow75)

## neutrophil percent by cbc_study_day and progressor status, with interaction
lmer.perc_neutrophils.cbc_study_day.progressor_split_cpep_rate_fast25_slow75.inc_interaction <-
  lmer(
    neutrophils ~
      cbc_study_day * progressor_split_cpep_rate_fast25_slow75 + (1|participant_id),
    data=master.cbc.merged[!duplicated(master.cbc.merged$cbc_visit_id_date),])
summary(lmer.perc_neutrophils.cbc_study_day.progressor_split_cpep_rate_fast25_slow75.inc_interaction)


##### plot % neutrophils over time by fast/slow progressor status #####

pdf("Fig_3A.pdf", w=10,h=6)
ggplot(
  data=
    master.cbc.merged[
      !duplicated(master.cbc.merged$cbc_visit_id_date) &
        !is.na(master.cbc.merged$progressor_split_cpep_rate_fast25_slow75) &
        !is.na(master.cbc.merged$neutrophils),] %>%
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
    data=master.cbc.merged[!duplicated(master.cbc.merged$cbc_visit_id_date),])
summary(lmer.perc_lymphocytes.cbc_study_day.progressor_split_cpep_rate_fast25_slow75)

## neutrophil percent by cbc_study_day and progressor status, with interaction
lmer.perc_lymphocytes.cbc_study_day.progressor_split_cpep_rate_fast25_slow75.inc_interaction <-
  lmer(
    lymphocytes ~
      cbc_study_day * progressor_split_cpep_rate_fast25_slow75 + (1|participant_id),
    data=master.cbc.merged[!duplicated(master.cbc.merged$cbc_visit_id_date),])
summary(lmer.perc_lymphocytes.cbc_study_day.progressor_split_cpep_rate_fast25_slow75.inc_interaction)


##### plot % lymphocytes over time by fast/slow progressor status #####

pdf("Fig_S9A.pdf", w=10, h=6)
ggplot(
  data=
    master.cbc.merged[
      !duplicated(master.cbc.merged$cbc_visit_id_date) &
        !is.na(master.cbc.merged$progressor_split_cpep_rate_fast25_slow75) &
        !is.na(master.cbc.merged$lymphocytes),] %>%
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
    data=master.cbc.merged[!duplicated(master.cbc.merged$cbc_visit_id_date),])
summary(lmer.lymphocytes_neutrophils_ratio.cbc_study_day.progressor_split_cpep_rate_fast25_slow75)

## lymphocytes_neutrophils_ratio by cbc_study_day and progressor status
lmer.lymphocytes_neutrophils_ratio.cbc_study_day.progressor_split_cpep_rate_fast25_slow75.inc_interaction <-
  lmer(
    lymphocytes_neutrophils_ratio ~
      cbc_study_day * progressor_split_cpep_rate_fast25_slow75 + (1|participant_id),
    data=master.cbc.merged[!duplicated(master.cbc.merged$cbc_visit_id_date),])
summary(lmer.lymphocytes_neutrophils_ratio.cbc_study_day.progressor_split_cpep_rate_fast25_slow75.inc_interaction)


##### plot lymphocyte-to-neutrophil ratio over time by fast/slow progressor status #####

data.tmp <-
  master.cbc.merged[
    !duplicated(master.cbc.merged$cbc_visit_id_date) &
      !is.na(master.cbc.merged$progressor_split_cpep_rate_fast25_slow75) &
      !is.na(master.cbc.merged$lymphocytes_neutrophils_ratio),] %>%
  cbind(
    lymphocytes_neutrophils_ratio_predicted =
      predict(lmer.lymphocytes_neutrophils_ratio.cbc_study_day.progressor_split_cpep_rate_fast25_slow75))

pdf("Fig_S9B.pdf", w=10, h=6)
ggplot(
  data=
    master.cbc.merged[
      !duplicated(master.cbc.merged$cbc_visit_id_date) &
        !is.na(master.cbc.merged$progressor_split_cpep_rate_fast25_slow75) &
        !is.na(master.cbc.merged$lymphocytes_neutrophils_ratio),] %>%
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


##### remove excess time columns from master.cbc.merged #####

master.cbc.merged <-
  master.cbc.merged[,!str_detect(colnames(master.cbc.merged), "scaled_")]


##### combine cell subset models into a single object for easier storage #####

cbc_models.all <- list()
for (i in ls_grep("^lm")) cbc_models.all[[i]] <- get(i)

rm_tmp(ask=FALSE)


##### filter and transform counts object and master.cbc.merged object #####

counts.filtered_cbcs <-
  calc_norm_counts(
    counts=counts.final, design=master.cbc.merged, libID_col="libid",
    min_cpm=0, min_libs_perc = 0,
    normalize=FALSE,
    return_DGEcounts=FALSE)
master.cbc.merged.filtered_cbcs <-
  master.cbc.merged[match(colnames(counts.filtered_cbcs), master.cbc.merged$libid),]


##### save objects for downstream use #####

save(file="T1D_placebos_data_2_for_downstream_analyses.RData",
     list=c(
       "master.cbc.merged", "master.cbc.merged.filtered_cbcs", "counts.filtered_cbcs",
       "counts.final", "auc2hr_models.all", "cbc_models.all", "vwts.all"))


