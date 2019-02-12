## scripts to run gene set analyses on T1D placebos project

##### set up environment: load packages #####

# load general packages
library(xlsx)
library(tidyverse)
library(colorspace)
library(RColorBrewer) 
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

# load analysis_specific packages
library(edgeR)
library(limma)
library(ordinal)
library(gdata) # needed for loading .xls data
library(lmerTest)

# load my packages of helper functions
library(geneSetTools)
library(miscHelpers)
library(countSubsetNorm)
library(limmaTools)


##### load/save data from previous scripts #####

load("T1D_placebos_data_1_for_downstream_analyses.RData")
load("T1D_placebos_data_2_for_downstream_analyses.RData")
load("T1D_placebos_data_3_for_downstream_analyses.RData")


##### calculate baseline median gene set expression for each library #####

median_gene_set_cols.tmp <-
  grep("^median", colnames(master.final), value=TRUE)
master.final[
  , paste0("baseline_", median_gene_set_cols.tmp)] <-
  as.numeric(NA)

for (i in unique(master.final$participant_id)) {
  baselines.tmp <-
    master.final[
      master.final$participant_id==i & master.final$rnaseq_baseline_visit,
      median_gene_set_cols.tmp]
  
  if (nrow(baselines.tmp) == 1)
    master.final[
      master.final$participant_id==i,
      paste0("baseline_", median_gene_set_cols.tmp)] <-
    baselines.tmp
}

rm_tmp(ask=FALSE)


##### calculate correlation of median CD19.mod and CXCR1.mod #####

master.final[
  master.final$rnaseq_baseline_visit,
  c("median_CD19.mod.removed_study_rna_batch", "median_CXCR1.mod.removed_study_rna_batch")] %>%
  cor(method="pearson", use="pairwise")


##### plot median expression at baseline of CXCR1.mod, CD19.mod, MPO.mod, and union of CXCR1.mod & MPO.mod #####

## plot each as a fit vs. age, separate plots, with a logarithmic function fit
pdf("Fig_4B.pdf", w=8, h=6)
ggplot(
  master.final[master.final$rnaseq_baseline_visit,],
  mapping=aes(x=age_years, y=median_CD19.mod.removed_study_rna_batch)) +
  geom_smooth(method="lm", formula=y~I(log(x)), color="black") +
  geom_point(size=2) +
  labs(x = "Age at diagnosis", y = "B cell gene set expression\n(median CD19.mod)")
dev.off()

pdf("Fig_3C.pdf", w=8, h=6)
ggplot(
  master.final[master.final$rnaseq_baseline_visit,],
  mapping=aes(x=age_years,  y=median_CXCR1.mod.removed_study_rna_batch)) +
  geom_smooth(method="lm", formula=y~I(log(x)), color="black") +
  geom_point(size=2) +
  labs(x = "Age at diagnosis", y = "Neutrophil gene set expression\n(median CXCR1.mod)")
dev.off()

pdf("Fig_3D.pdf", w=8, h=6)
ggplot(
  master.final[master.final$rnaseq_baseline_visit,],
  mapping=aes(x=age_years,  y=median_MPO.mod.removed_study_rna_batch)) +
  geom_smooth(method="lm", formula=y~I(log(x)), color="black") +
  geom_point(size=2) +
  labs(x = "Age at diagnosis", y = "Neutrophil gene set expression\n(median MPO.mod)")
dev.off()


##### fit models to median gene set expression vs. log_age_years #####

## fit model to median CD19.mod.removed_study_rna_batch
lm.median_CD19.mod.removed_study_rna_batch_vs_log_age_years <-
  lm(median_CD19.mod.removed_study_rna_batch ~ log(age_years),
     data=master.final[master.final$rnaseq_baseline_visit,])
summary(lm.median_CD19.mod.removed_study_rna_batch_vs_log_age_years) # p=4.5E-8

## fit model to median CXCR1.mod.removed_study_rna_batch
lm.median_CXCR1.mod.removed_study_rna_batch_vs_log_age_years <-
  lm(median_CXCR1.mod.removed_study_rna_batch ~ log(age_years),
     data=master.final[master.final$rnaseq_baseline_visit,])
summary(lm.median_CXCR1.mod.removed_study_rna_batch_vs_log_age_years) # p=5.0E-3

## fit model to median MPO.mod.removed_study_rna_batch
lm.median_MPO.mod.removed_study_rna_batch_vs_log_age_years <-
  lm(median_MPO.mod.removed_study_rna_batch ~ log(age_years),
     data=master.final[master.final$rnaseq_baseline_visit,])
summary(lm.median_MPO.mod.removed_study_rna_batch_vs_log_age_years) # p=0.19


##### plot and check correlations between neutrophils/lymphocytes from CBCs, and expression of gene modules #####

## plots of neutrophils from CBCs vs. neutrophil module gene expression

# for CXCR1.mod
pdf("Fig_S7A.pdf", w=7, h=6)
ggplot(
  master.cbc.merged,
  mapping=aes(x=median_CXCR1.mod.removed_study_rna_batch, y=neutrophils)) +
  geom_point(size=2)+
  geom_smooth(method="lm", linetype="dashed", color="black", se=FALSE) +
  labs(x="Neutrophil gene set expression\n(median CXCR1.mod)",
       y="% neutrophils by CBC")
dev.off()

# for MPO.mod
pdf("Fig_S7A.pdf", w=7, h=6)
ggplot(
  master.cbc.merged,
  mapping=aes(x=median_MPO.mod.removed_study_rna_batch, y=neutrophils)) +
  geom_point(size=2)+
  geom_smooth(method="lm", linetype="dashed", color="black", se=FALSE) +
  labs(x="Neutrophil gene set expression\n(median MPO.mod)",
       y="% neutrophils by CBC")
dev.off()

# correlations between neutrophil percent from CBC, and neutrophil module expression from RNAseq
cor.test(master.cbc.merged$neutrophils,
         master.cbc.merged$median_CXCR1.mod.removed_study_rna_batch,
         method="pearson")

cor.test(master.cbc.merged$neutrophils,
         master.cbc.merged$median_MPO.mod.removed_study_rna_batch,
         method="pearson")


##### fit models to median gene set expression by progressor status and rnaseq_study_day #####

## for CD19.mod
lmer.median_CD19.mod.removed_study_rna_batch.rnaseq_study_day.progressor_split_cpep_rate_fast25_slow75 <-
  lmer(
    median_CD19.mod.removed_study_rna_batch ~
      rnaseq_study_day + progressor_split_cpep_rate_fast25_slow75 + (1|participant_id),
    data=master.final,
    na.action=na.omit)
summary(lmer.median_CD19.mod.removed_study_rna_batch.rnaseq_study_day.progressor_split_cpep_rate_fast25_slow75)

## for CD19.mod, with interaction
lmer.median_CD19.mod.removed_study_rna_batch.rnaseq_study_day.progressor_split_cpep_rate_fast25_slow75.inc_interaction <-
  lmer(
    median_CD19.mod.removed_study_rna_batch ~
      rnaseq_study_day * progressor_split_cpep_rate_fast25_slow75 + (1|participant_id),
    data=master.final)
summary(lmer.median_CD19.mod.removed_study_rna_batch.rnaseq_study_day.progressor_split_cpep_rate_fast25_slow75.inc_interaction)

## for CXCR1.mod
lmer.median_CXCR1.mod.removed_study_rna_batch.rnaseq_study_day.progressor_split_cpep_rate_fast25_slow75 <-
  lmer(
    median_CXCR1.mod.removed_study_rna_batch ~
      rnaseq_study_day + progressor_split_cpep_rate_fast25_slow75 + (1|participant_id),
    data=master.final,
    na.action=na.omit)
summary(lmer.median_CXCR1.mod.removed_study_rna_batch.rnaseq_study_day.progressor_split_cpep_rate_fast25_slow75)

## for CXCR1.mod, with interaction
lmer.median_CXCR1.mod.removed_study_rna_batch.rnaseq_study_day.progressor_split_cpep_rate_fast25_slow75.inc_interaction <-
  lmer(
    median_CXCR1.mod.removed_study_rna_batch ~
      rnaseq_study_day * progressor_split_cpep_rate_fast25_slow75 + (1|participant_id),
    data=master.final)
summary(lmer.median_CXCR1.mod.removed_study_rna_batch.rnaseq_study_day.progressor_split_cpep_rate_fast25_slow75.inc_interaction)

## for MPO.mod
lmer.median_MPO.mod.removed_study_rna_batch.rnaseq_study_day.progressor_split_cpep_rate_fast25_slow75 <-
  lmer(
    median_MPO.mod.removed_study_rna_batch ~
      rnaseq_study_day + progressor_split_cpep_rate_fast25_slow75 + (1|participant_id),
    data=master.final,
    na.action=na.omit)
summary(lmer.median_MPO.mod.removed_study_rna_batch.rnaseq_study_day.progressor_split_cpep_rate_fast25_slow75)

## for MPO.mod, with interaction
lmer.median_MPO.mod.removed_study_rna_batch.rnaseq_study_day.progressor_split_cpep_rate_fast25_slow75.inc_interaction <-
  lmer(
    median_MPO.mod.removed_study_rna_batch ~
      rnaseq_study_day * progressor_split_cpep_rate_fast25_slow75 + (1|participant_id),
    data=master.final)
summary(lmer.median_MPO.mod.removed_study_rna_batch.rnaseq_study_day.progressor_split_cpep_rate_fast25_slow75.inc_interaction)


##### plot median gene set expression over time, by progressor status #####

# plot median gene expression over time by fast/slow progressor status

gene_set_summaries.tmp <-
  c("median_CD19.mod.removed_study_rna_batch",
    "median_CXCR1.mod.removed_study_rna_batch",
    "median_MPO.mod.removed_study_rna_batch")
axis_labs.tmp <-
  c("median_CD19.mod.removed_study_rna_batch" = "B cell gene set expression\n(median CD19.mod)",
    "median_CXCR1.mod.removed_study_rna_batch" = "Neutrophil gene set expression\n(median CXCR1.mod)",
    "median_MPO.mod.removed_study_rna_batch" = "Neutrophil gene set expression\n(median MPO.mod)")
filename.tmp <-
  c("median_CD19.mod.removed_study_rna_batch" = "Fig_4A.pdf",
    "median_CXCR1.mod.removed_study_rna_batch" = "Fig_S6A.pdf",
    "median_MPO.mod.removed_study_rna_batch" = "Fig_S6B.pdf")

master.tmp <- master.final[!is.na(master.final$rnaseq_study_day),]

for (i in gene_set_summaries.tmp) {
  lm.tmp <-
    lm(
      formula(paste(i, "~ rnaseq_study_day + progressor_split_cpep_rate_fast25_slow75")),
      data=master.tmp)
  pdf(filename.tmp[i], w=10, h=6)
  print(
    ggplot() +
      # geom_jitter(width=2, height=0) +
      geom_point(
        master.tmp,
        mapping=aes_string(
          x="rnaseq_study_day", y=i,
          group="progressor_split_cpep_rate_fast25_slow75",
          color="progressor_split_cpep_rate_fast25_slow75"),
        alpha=0.7, size=3, shape=16) +
      geom_smooth(
        mapping=aes(
          x=master.tmp$rnaseq_study_day[
            master.tmp$progressor_split_cpep_rate_fast25_slow75 %in% "fast"],
          y=predict(lm.tmp)[
            master.tmp$progressor_split_cpep_rate_fast25_slow75 %in% "fast"]),
        method="lm", span=0.7, size=2.0, se=FALSE, color="black") +
      geom_smooth(
        mapping=aes(
          x=master.tmp$rnaseq_study_day[
            master.tmp$progressor_split_cpep_rate_fast25_slow75 %in% "fast"],
          y=predict(lm.tmp)[
            master.tmp$progressor_split_cpep_rate_fast25_slow75 %in% "fast"]),
        method="lm", size=1.5, se=FALSE, color="red") +
      geom_smooth(
        mapping=aes(
          x=master.tmp$rnaseq_study_day[
            master.tmp$progressor_split_cpep_rate_fast25_slow75 %in% "slow"],
          y=predict(lm.tmp)[
            master.tmp$progressor_split_cpep_rate_fast25_slow75 %in% "slow"]),
        method="lm", span=0.7, size=2.0, se=FALSE, color="black") +
      geom_smooth(
        mapping=aes(
          x=master.tmp$rnaseq_study_day[
            master.tmp$progressor_split_cpep_rate_fast25_slow75 %in% "slow"],
          y=predict(lm.tmp)[
            master.tmp$progressor_split_cpep_rate_fast25_slow75 %in% "slow"]),
        method="lm", size=1.5, se=FALSE, color="blue") +
      scale_color_manual(name="progressor", values=c("red", "blue")) +
      scale_x_continuous(breaks=seq(0,1500, by=250)) +
      labs(x="Days in study",
           y=axis_labs.tmp[i]))
  dev.off()
}

rm_tmp(ask=FALSE)


##### fit model to median gene set expression vs. cpeptide_slope, by age tertiles #####

## for CD19.mod
lm.cpeptide_slope_vs_median_CD19.mod.removed_study_rna_batch.age_split3.with_interaction.t0 <-
  lm(log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept ~
       median_CD19.mod.removed_study_rna_batch +
       age_split3 +
       median_CD19.mod.removed_study_rna_batch * age_split3,
     data=
       master.final[
         with(master.final,
              !is.na(log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept) &
                rnaseq_baseline_visit),])
summary(lm.cpeptide_slope_vs_median_CD19.mod.removed_study_rna_batch.age_split3.with_interaction.t0)

## for CXCR1.mod
lm.cpeptide_slope_vs_median_CXCR1.mod.removed_study_rna_batch.age_split3.with_interaction.t0 <-
  lm(log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept ~
       median_CXCR1.mod.removed_study_rna_batch +
       age_split3 +
       median_CXCR1.mod.removed_study_rna_batch * age_split3,
     data=
       master.final[
         with(master.final,
              !is.na(log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept) &
                rnaseq_baseline_visit),])
summary(lm.cpeptide_slope_vs_median_CXCR1.mod.removed_study_rna_batch.age_split3.with_interaction.t0)

## for MPO.mod
lm.cpeptide_slope_vs_median_MPO.mod.removed_study_rna_batch.age_split3.with_interaction.t0 <-
  lm(log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept ~
       median_MPO.mod.removed_study_rna_batch +
       age_split3 +
       median_MPO.mod.removed_study_rna_batch * age_split3,
     data=
       master.final[
         with(master.final,
              !is.na(log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept) &
                rnaseq_baseline_visit),])
summary(lm.cpeptide_slope_vs_median_MPO.mod.removed_study_rna_batch.age_split3.with_interaction.t0)


##### plot C-peptide rate vs. median gene set expression, by age tertile #####

## for CXCR1.mod
pdf("Fig_S10.pdf", w=12, h=7)
ggplot(
  master.final %>%
    dplyr::filter(rnaseq_baseline_visit),
  mapping=aes(
    x=median_CXCR1.mod.removed_study_rna_batch,
    y=log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept,
    color=age_split3)) +
  geom_smooth(
    method="lm", se=FALSE) +
  geom_point(size=2) +
  labs(x="neutrophil gene expression\n(median CXCR1.mod)",
       y="Rate of C-peptide change") +
  scale_color_manual(
    "Age group", values=ggthemes::colorblind_pal()(4)[-1],
    labels=c(
      "Youngest third (8.5 to 13.4)",
      "Middle third (13.4 to 17.3)",
      "Oldest third (17.3 to 46.2)"))
dev.off()

## for CD19.mod
pdf("Fig_4C.pdf", w=12, h=7)
ggplot(
  master.final %>%
    dplyr::filter(rnaseq_baseline_visit),
  mapping=aes(
    x=median_CD19.mod.removed_study_rna_batch,
    y=log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept,
    color=age_split3)) +
  geom_smooth(
    method="lm", se=FALSE) +
  geom_point(size=2) +
  labs(x="B cell gene expression\n(median CD19.mod)",
       y="Rate of C-peptide change") +
  scale_color_manual(
    "Age group", values=ggthemes::colorblind_pal()(4)[-1],
    labels=c(
      "Youngest third (8.5 to 13.4)",
      "Middle third (13.4 to 17.3)",
      "Oldest third (17.3 to 46.2)"))
dev.off()


##### sliding age window of median CD19.mod gene expression vs. rate of C-peptide change #####

# filter to only baseline visits (also cuts to 1 visit per subject), and arrange by age
master.t0 <-
  master.final %>%
  dplyr::filter(rnaseq_baseline_visit) %>%
  arrange(age_years)

# set window size
window_size.CD19.mod_cpeptide_slope <- 35

cor.pearson.sliding_window.median_CD19.mod.removed_study_rna_batch.cpeptide_slope.t0 <-
  data.frame(
    window_start=1:(nrow(master.t0)-window_size.CD19.mod_cpeptide_slope+1),
    age_min=as.numeric(NA),
    age_max=as.numeric(NA),
    median_age=as.numeric(NA),
    correlation=as.numeric(NA),
    p_value=as.numeric(NA))
for (i in seq.int(nrow(cor.pearson.sliding_window.median_CD19.mod.removed_study_rna_batch.cpeptide_slope.t0))) {
  start.tmp <-
    cor.pearson.sliding_window.median_CD19.mod.removed_study_rna_batch.cpeptide_slope.t0$window_start[i]
  master.tmp <-
    master.t0[start.tmp:(start.tmp+window_size.CD19.mod_cpeptide_slope-1),]
  cor.pearson.sliding_window.median_CD19.mod.removed_study_rna_batch.cpeptide_slope.t0[
    i,c("age_min", "age_max")] <-
    range(master.tmp$age_years, na.rm=TRUE)
  cor.pearson.sliding_window.median_CD19.mod.removed_study_rna_batch.cpeptide_slope.t0$median_age[i] <-
    median(master.tmp$age_years, na.rm=TRUE)
  cor.pearson.sliding_window.median_CD19.mod.removed_study_rna_batch.cpeptide_slope.t0$correlation[i] <-
    cor(
      master.tmp$median_CD19.mod.removed_study_rna_batch,
      master.tmp$log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept,
      method="pearson", use="pairwise")
  cor.pearson.sliding_window.median_CD19.mod.removed_study_rna_batch.cpeptide_slope.t0$p_value[i] <-
    cor.test(
      master.tmp$median_CD19.mod.removed_study_rna_batch,
      master.tmp$log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept,
      method="pearson")[["p.value"]]
}

# now plot the results
pdf("Fig_S11.pdf", w=8, h=7)
ggplot(
  cor.pearson.sliding_window.median_CD19.mod.removed_study_rna_batch.cpeptide_slope.t0,
  mapping=aes(
    x=median_age,
    y=correlation)) +
  geom_line() +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(x="Median age of window")
dev.off()

rm_tmp(ask=FALSE)


##### plot C-peptide rate vs. median B cell gene set expression, by age split at sliding window age break #####

age_split_sliding_window <- 20
master.tmp <-
  master.final %>%
  dplyr::filter(rnaseq_baseline_visit) %>%
  mutate(age_split_sliding_window =
           factor(ifelse(age_years < age_split_sliding_window, "younger", "older"),
                  levels=c("younger", "older")))
pdf("Fig_4D.pdf", w=11.5, h=7)
ggplot(
  master.tmp,
  mapping=aes(
    x=median_CD19.mod.removed_study_rna_batch,
    y=log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept,
    color=age_split_sliding_window)) +
  geom_smooth(
    method="lm", se=FALSE) +
  geom_point(size=2.5) +
  labs(x="B cell gene expression\n(median CD19.mod)",
       y="Rate of C-peptide change") +
  scale_color_manual(
    "Age group", values=ggthemes::colorblind_pal()(4)[c(2,4)],
    labels=c(
      "Younger (< 20 years)",
      "Older (> 20 years)"))
dev.off()

lm(log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept ~
     age_split_sliding_window * median_CD19.mod.removed_study_rna_batch,
   data=master.tmp) %>%
  summary()

rm_tmp(ask=FALSE)


##### save objects for downstream use #####

save(file="T1D_placebos_data_4_for_downstream_analyses.RData",
     list=c("master.final", "master.t0",
            ls_grep("^lm")))
