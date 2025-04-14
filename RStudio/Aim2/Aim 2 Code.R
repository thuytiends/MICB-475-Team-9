##### MICB 475 HIV Project Aim 2 - Core microbiome, indicator analysis, Deseq 2

library(tidyverse)
library(DESeq2)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(indicspecies)
library(dplyr)

#### Load phyloseq object in and create substratifed phyloseq objects for downstream analysis

load("hiv_rare.RData")
hiv_rare <- subset_samples(hiv_rare,Visit_Cat == "2nd Visit")

labor_rename <- as.data.frame(sample_data(hiv_rare))
labor_rename$Visit_Age <- ifelse(labor_rename$Visit_Age > 45, "Older than 45", "Younger than 45")
labor_rename$Labor_category <- ifelse(labor_rename$Labor_category %in% c("Non_Labor_Houseworker","Non_Labor_Other","Unemployed"), "Non-Labor", "Labor")
sample_data(hiv_rare) <- sample_data(labor_rename)

hiv_positive <- subset_samples(hiv_rare, HIV_Status == "Positive")
hiv_negative <- subset_samples(hiv_rare, HIV_Status == "Negative")

#### CORE MICROBIOME ANALYSIS

## Convert to relative abundance
hiv_positive_RA <- transform_sample_counts(hiv_positive, fun=function(x) x/sum(x))
hiv_negative_RA <- transform_sample_counts(hiv_negative, fun=function(x) x/sum(x))

hiv_positive_prim_edu <- subset_samples(hiv_positive_RA, `Education_Level`=="Primary")
hiv_positive_second_edu <- subset_samples(hiv_positive_RA, `Education_Level`=="Secondary")
hiv_positive_tert_edu <- subset_samples(hiv_positive_RA, `Education_Level`=="Tertiary")

hiv_positive_old <- subset_samples(hiv_positive_RA, `Visit_Age`== "Older than 45") 
hiv_positive_young <- subset_samples(hiv_positive_RA, `Visit_Age`== "Younger than 45") 

hiv_positive_lab <- subset_samples(hiv_positive_RA, `Labor_category`== "Labor")
hiv_positive_non_lab <- subset_samples(hiv_positive_RA, `Labor_category`== "Non-Labor")

hiv_negative_prim_edu <- subset_samples(hiv_negative_RA, `Education_Level`=="Primary")
hiv_negative_second_edu <- subset_samples(hiv_negative_RA, `Education_Level`=="Secondary")
hiv_negative_tert_edu <- subset_samples(hiv_negative_RA, `Education_Level`=="Tertiary")

hiv_negative_old <- subset_samples(hiv_negative_RA, `Visit_Age`== "Older than 45") 
hiv_negative_young <- subset_samples(hiv_negative_RA, `Visit_Age`== "Younger than 45") 

hiv_negative_lab <- subset_samples(hiv_negative_RA, `Labor_category`== "Labor")
hiv_negative_non_lab <- subset_samples(hiv_negative_RA, `Labor_category`== "Non-Labor")

## What ASVs are found in more than 70% of samples in each category? Trying changing the prevalence to see what happens

hiv_positive_prim_edu_ASVs <- core_members(hiv_positive_prim_edu, detection=0, prevalence = 0.7)
hiv_positive_sec_edu_ASVs <- core_members(hiv_positive_second_edu, detection=0, prevalence = 0.7)
hiv_positive_tert_edu_ASVs <- core_members(hiv_positive_tert_edu, detection=0, prevalence = 0.7)

hiv_positive_young_ASVs <- core_members(hiv_positive_young, detection=0, prevalence = 0.7)
hiv_positive_old_ASVs <- core_members(hiv_positive_old, detection=0, prevalence = 0.7)

hiv_positive_lab_ASVs <- core_members(hiv_positive_lab, detection=0, prevalence = 0.7)
hiv_positive_non_lab_ASVs <- core_members(hiv_positive_non_lab, detection=0, prevalence = 0.7)

hiv_negative_prim_edu_ASVs <- core_members(hiv_negative_prim_edu, detection=0, prevalence = 0.7)
hiv_negative_sec_edu_ASVs <- core_members(hiv_negative_second_edu, detection=0, prevalence = 0.7)
hiv_negative_tert_edu_ASVs <- core_members(hiv_negative_tert_edu, detection=0, prevalence = 0.7)

hiv_negative_young_ASVs <- core_members(hiv_negative_young, detection=0, prevalence = 0.7)
hiv_negative_old_ASVs <- core_members(hiv_negative_old, detection=0, prevalence = 0.7)

hiv_negative_lab_ASVs <- core_members(hiv_negative_lab, detection=0, prevalence = 0.7)
hiv_negative_non_lab_ASVs <- core_members(hiv_negative_non_lab, detection=0, prevalence = 0.7)

## What are these ASVs?

tax_table(prune_taxa(hiv_positive_prim_edu_ASVs,hiv_rare))
tax_table(prune_taxa(hiv_positive_sec_edu_ASVs,hiv_rare))
tax_table(prune_taxa(hiv_positive_tert_edu_ASVs,hiv_rare))
tax_table(prune_taxa(hiv_positive_young_ASVs,hiv_rare))
tax_table(prune_taxa(hiv_positive_old_ASVs,hiv_rare))
tax_table(prune_taxa(hiv_positive_lab_ASVs,hiv_rare))
tax_table(prune_taxa(hiv_positive_non_lab_ASVs,hiv_rare))

tax_table(prune_taxa(hiv_negative_prim_edu_ASVs,hiv_rare))
tax_table(prune_taxa(hiv_negative_sec_edu_ASVs,hiv_rare))
tax_table(prune_taxa(hiv_negative_tert_edu_ASVs,hiv_rare))
tax_table(prune_taxa(hiv_negative_young_ASVs,hiv_rare))
tax_table(prune_taxa(hiv_negative_old_ASVs,hiv_rare))
tax_table(prune_taxa(hiv_negative_lab_ASVs,hiv_rare))
tax_table(prune_taxa(hiv_negative_non_lab_ASVs,hiv_rare))

#HIV Positive

##Can plot those ASVs' relative abundance 
## Is there a way to average the ASVs?
prune_taxa(hiv_positive_prim_edu_ASVs,hiv_positive_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Education_Level`, scales ="free")

prune_taxa(hiv_positive_sec_edu_ASVs,hiv_positive_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Education_Level`, scales ="free")

prune_taxa(hiv_positive_tert_edu_ASVs,hiv_positive_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Education_Level`, scales ="free")

prune_taxa(hiv_positive_young_ASVs,hiv_positive_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Visit_Age`, scales ="free")

prune_taxa(hiv_positive_old_ASVs,hiv_positive_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Visit_Age`, scales ="free")

prune_taxa(hiv_positive_lab_ASVs,hiv_positive_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Labor_category`, scales ="free")

prune_taxa(hiv_positive_non_lab_ASVs,hiv_positive_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Labor_category`, scales ="free")

## What if we wanted to make a Venn diagram of all the ASVs that showed up in each treatment?
hiv_positive_lab_list <- core_members(hiv_positive_lab, detection=0.001, prevalence = 0.10)
hiv_positive_non_lab_list <- core_members(hiv_positive_non_lab, detection=0.001, prevalence = 0.10)

hiv_negative_lab_list <- core_members(hiv_negative_lab, detection=0.001, prevalence = 0.10)
hiv_negative_non_lab_list <- core_members(hiv_negative_non_lab, detection=0.001, prevalence = 0.10)

hiv_positive_old_list <- core_members(hiv_positive_old, detection=0.001, prevalence = 0.10)
hiv_positive_young_list <- core_members(hiv_positive_young, detection=0.001, prevalence = 0.10)

hiv_positive_prim_edu_list <- core_members(hiv_positive_prim_edu, detection=0.001, prevalence = 0.10)
hiv_positive_sec_edu_list <- core_members(hiv_positive_second_edu, detection=0.001, prevalence = 0.10)
hiv_positive_tert_edu_list <- core_members(hiv_positive_tert_edu, detection=0.001, prevalence = 0.10)



hiv_neg_pos_lab_list_full <- list(HIV_Positive = hiv_positive_lab_list, HIV_Negative = hiv_negative_lab_list)
hiv_lab_nolab_list_full <- list(Labour = hiv_positive_lab_list, No_Labour = hiv_positive_non_lab_list)
hiv_old_young_list_full <- list(Older_than_45 = hiv_positive_old_list, Younger_than_45 = hiv_positive_young_list)
hiv_education_full <- list(Primary = hiv_positive_prim_edu_list, Secondary = hiv_positive_sec_edu_list, Third = hiv_positive_tert_edu_list)


library("sf")

# Create a Venn diagram using all the ASVs 
hiv_neg_pos_lab_venn <- ggVennDiagram(x = hiv_neg_pos_lab_list_full)

hiv_pos_lab_nonlab_venn <- ggVennDiagram(x = hiv_lab_nolab_list_full)
hiv_pos_young_old_venn <- ggVennDiagram(x = hiv_old_young_list_full)
#younger on top older on bottom
hiv_pos_education_venn <- ggVennDiagram(x = hiv_education_full) 

ggsave("hiv_neg_pos_labour_venn.png", hiv_neg_pos_lab_venn)
ggsave("hiv_pos_labour_nonlabour_venn.png", hiv_pos_lab_nonlab_venn)
ggsave("hiv_pos_young_aging_venn.png", hiv_pos_young_old_venn)
ggsave("hiv_pos_education_venn.png", hiv_pos_education_venn)

#Venn diagram labels are off, need to figure out a way to recenter them and refine the diagrams

# ------------------------------------------------------------------------------

#### INDICATOR SPECIES ANALYSIS

# glom to Genus HIV Positive

hiv_positive_genus <- tax_glom(hiv_positive, "Genus", NArm = FALSE)
hiv_positive_genus_RA <- transform_sample_counts(hiv_positive_genus, fun=function(x) x/sum(x))

###HIV POSITIVE###

#Education Level
isa_hiv_positive_edu <- multipatt(t(otu_table(hiv_positive_genus_RA)), cluster = sample_data(hiv_positive_genus_RA)$`Education_Level`)

summary(isa_hiv_positive_edu)

taxtable <- tax_table(hiv_rare) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# consider that your table is only going to be resolved up to the genus level, be wary of anything beyond the glomed taxa level stat closer to 1, better indicator

isa_hiv_positive_edu$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% View()

#Got 3 positive ASVs that were significant that were in tertiary

#Labour

isa_hiv_positive_labour <- multipatt(t(otu_table(hiv_positive_genus_RA)), cluster = sample_data(hiv_positive_genus_RA)$`Labor_category`)

summary(isa_hiv_positive_labour)

taxtable <- tax_table(hiv_rare) %>% as.data.frame() %>% rownames_to_column(var="ASV")

isa_hiv_positive_labour$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% View()

#Age

isa_hiv_positive_age <- multipatt(t(otu_table(hiv_positive_genus_RA)), cluster = sample_data(hiv_positive_genus_RA)$`Visit_Age`)

summary(isa_hiv_positive_age)

taxtable <- tax_table(hiv_rare) %>% as.data.frame() %>% rownames_to_column(var="ASV")

isa_hiv_positive_age$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% View()

# ------------------------------------------------------------------------------

#### DESEQ-2 ANALYSIS
## DESEQ, NOTE: If you get a zeros error, then you need to add '1' count to all reads

### Load data ###
load("hiv_final.RData")

#metadata reorganizing
hiv_final <- subset_samples(hiv_final,Visit_Cat == "2nd Visit")

labor_rename <- as.data.frame(sample_data(hiv_final))
labor_rename$Visit_Age <- ifelse(labor_rename$Visit_Age > 45, "Older than 45", "Younger than 45")
labor_rename$Labor_category <- ifelse(labor_rename$Labor_category %in% c("Non_Labor_Houseworker","Non_Labor_Other","Unemployed"), "Non-Labor", "Labor")
sample_data(hiv_final) <- sample_data(labor_rename)

hiv_pos <- subset_samples(hiv_final, HIV_Status == "Positive")
hiv_neg <- subset_samples(hiv_final, HIV_Status == "Negative")
head(sample_data(hiv_final))

colnames(sample_data(hiv_final))

#Age
hiv_deseq_agegrp <- phyloseq_to_deseq2(hiv_final, ~`Visit_Age`)
hiv_agegrp_plus1 <- transform_sample_counts(hiv_pos, function(x) x+1)

hiv_deseq_agegrp <- phyloseq_to_deseq2(hiv_agegrp_plus1, ~`Visit_Age`)
DESEQ_hiv_agegrp <- DESeq(hiv_deseq_agegrp)

results_hiv_agegrp <- results(DESEQ_hiv_agegrp, tidy=TRUE, 
                              #this will ensure that No is your reference group
                              contrast = c("Visit_Age","Older than 45","Younger than 45"))

View(results_hiv_agegrp)
colnames(colData(DESEQ_hiv_labor))
colnames(colData(DESEQ_hiv_agegrp))
#Education
hiv_deseq_Education <- phyloseq_to_deseq2(hiv_pos, ~`Education_Level`)
hiv_Education_plus1 <- transform_sample_counts(hiv_pos, function(x) x+1)

hiv_deseq_Education <- phyloseq_to_deseq2(hiv_Education_plus1, ~`Education_Level`)
DESEQ_hiv_Education <- DESeq(hiv_deseq_Education)


results_hiv_Education <- results(DESEQ_hiv_Education, tidy=TRUE, 
                              #this will ensure that No is your reference group
                              contrast = c("Education_Level","Secondary","Tertiary"))

View(results_hiv_Education)

#Labor
hiv_deseq_labor <- phyloseq_to_deseq2(hiv_pos, ~`Labor_category`)
hiv_labor_plus1 <- transform_sample_counts(hiv_pos, function(x) x+1)

hiv_deseq_labor <- phyloseq_to_deseq2(hiv_labor_plus1, ~`Labor_category`)
DESEQ_hiv_labor <- DESeq(hiv_deseq_labor)


results_hiv_labor <- results(DESEQ_hiv_labor, tidy=TRUE, 
                                 #this will ensure that No is your reference group
                                 contrast = c("Labor_category","Labor","Non-Labor"))
levels(DESEQ_hiv_labor$Labor_category)

View(results_hiv_labor)
## Volcano plot: effect size VS significance
## Make variable to color by whether it is significant + large change
#Age
vol_plot_hiv_agegrp <- results_hiv_agegrp %>%   
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 2) %>%   
  ggplot() +   
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = significant)) +   
  scale_color_manual(values = c("FALSE" = "#fdae61", "TRUE" = "#d7191c")) +   
  theme_minimal() +  
  theme(
    axis.text = element_text(size = 16),      # Set axis text size
    axis.title = element_text(size = 18)      # Set axis title size
  )

vol_plot_hiv_agegrp

ggsave(filename="vol_plot_hiv_agegrp_hivpos.png",vol_plot_hiv_agegrp)

#Labor
vol_plot_hiv_labor <- results_hiv_labor %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 2) %>%
  ggplot() +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  scale_color_manual(values = c("FALSE" = "#fdae61", "TRUE" = "#d7191c")) +
  theme_minimal() +  
  theme(
    axis.text = element_text(size = 16),      # Set axis text size
    axis.title = element_text(size = 18)      # Set axis title size
  )

vol_plot_hiv_labor

ggsave(filename="vol_plot_hiv_labor_hivpos.png",vol_plot_hiv_labor)

#Education volcano plot
vol_plot_hiv_Education <- results_hiv_Education %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 2) %>%
  ggplot() +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  scale_color_manual(values = c("FALSE" = "#fdae61", "TRUE" = "#d7191c")) +
  theme_minimal() +  
  theme(
    axis.text = element_text(size = 16),      # Set axis text size
    axis.title = element_text(size = 18)      # Set axis title size
  )

vol_plot_hiv_Education

ggsave(filename="vol_plot_hiv_Education_sec_ter_hivpos.png",vol_plot_hiv_Education)


# To get table of results
sigASVs <- results_hiv_labor %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

sigASVs <- results_hiv_agegrp %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

sigASVs <- results_hiv_Education %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)


View(sigASVs)
nrow(sigASVs)
# Get only asv names
sigASVs_vec <- sigASVs %>%
  pull(ASV)

## Prune phyloseq file

#Age
hiv_agegrp_DESeq <- prune_taxa(sigASVs_vec,hiv_pos)
hiv_agegrp_sigASVs <- tax_table(hiv_agegrp_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

DESeqsig_hiv_agegrp <-ggplot(hiv_agegrp_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

DESeqsig_hiv_agegrp
ggsave(filename = "DESeqsig_hiv_agegrp_hivpos.png", plot = DESeqsig_hiv_agegrp, width = 12, height = 8, dpi = 300)

#Labor
hiv_labor_DESeq <- prune_taxa(sigASVs_vec,hiv_pos)
hiv_labor_sigASVs <- tax_table(hiv_labor_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

DESeqsig_hiv_labor <-ggplot(hiv_labor_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

DESeqsig_hiv_labor
ggsave(filename = "DESeqsig_hiv_labor_hivpos.png", plot = DESeqsig_hiv_labor, width = 12, height = 8, dpi = 300)


#Education
hiv_Education_DESeq <- prune_taxa(sigASVs_vec,hiv_pos)
hiv_Education_sigASVs <- tax_table(hiv_Education_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

DESeqsig_hiv_Education <-ggplot(hiv_Education_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

DESeqsig_hiv_Education
ggsave(filename = "DESeqsig_hiv_Education_pri_sec_hivpos.png", plot = DESeqsig_hiv_Education, width = 12, height = 8, dpi = 300)

