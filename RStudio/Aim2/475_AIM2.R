#### Load data ####
load("hiv_final.RData")

#metadata reorganizing
hiv_final <- subset_samples(hiv_final,Visit_Cat == "2nd Visit")

labor_rename <- as.data.frame(sample_data(hiv_final))
labor_rename$Visit_Age <- ifelse(labor_rename$Visit_Age > 45, "Older than 45", "Younger than 45")
labor_rename$Labor_category <- ifelse(labor_rename$Labor_category %in% c("Non_Labor_Houseworker","Non_Labor_Other","Unemployed"), "Non-Labor", "Labor")
sample_data(hiv_final) <- sample_data(labor_rename)

hiv_pos <- subset_samples(hiv_final, HIV_Status == "Positive")
hiv_neg <- subset_samples(hiv_final, HIV_Status == "Negative")

#################################################################################################
#### CORE MICROBIOME ####

library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

# Convert to relative abundance
hiv_pos_RA <- transform_sample_counts(hiv_pos, fun=function(x) x/sum(x))
hiv_neg_RA <- transform_sample_counts(hiv_neg, fun=function(x) x/sum(x))

hiv_pos_1edu <- subset_samples(hiv_pos_RA, `Education_Level`=="Primary")
hiv_pos_2edu <- subset_samples(hiv_pos_RA, `Education_Level`=="Secondary")
hiv_pos_3edu <- subset_samples(hiv_pos_RA, `Education_Level`=="Tertiary")

hiv_pos_aging <- subset_samples(hiv_pos_RA, `Visit_Age`== "Older than 45") ##18-45
hiv_pos_young <- subset_samples(hiv_pos_RA, `Visit_Age`== "Younger than 45") ##46-70

hiv_pos_lab <- subset_samples(hiv_pos_RA, `Labor_category`== "Labor")
hiv_pos_nonlab <- subset_samples(hiv_pos_RA, `Labor_category`== "Non-Labor")

hiv_neg_1edu <- subset_samples(hiv_neg_RA, `Education_Level`=="Primary")
hiv_neg_2edu <- subset_samples(hiv_neg_RA, `Education_Level`=="Secondary")
hiv_neg_3edu <- subset_samples(hiv_neg_RA, `Education_Level`=="Tertiary")

hiv_neg_aging <- subset_samples(hiv_neg_RA, `Visit_Age`== "Older than 45") ##18-45
hiv_neg_young <- subset_samples(hiv_neg_RA, `Visit_Age`== "Younger than 45") ##46-70

hiv_neg_lab <- subset_samples(hiv_neg_RA, `Labor_category`== "Labor")
hiv_neg_nonlab <- subset_samples(hiv_neg_RA, `Labor_category`== "Non-Labor")


# What ASVs are found in more than 70% of samples in each category?
# trying changing the prevalence to see what happens
hiv_pos_1edu_ASVs <- core_members(hiv_pos_1edu, detection=0, prevalence = 0.7)
hiv_pos_2edu_ASVs <- core_members(hiv_pos_2edu, detection=0, prevalence = 0.7)
hiv_pos_3edu_ASVs <- core_members(hiv_pos_3edu, detection=0, prevalence = 0.7)

hiv_pos_young_ASVs <- core_members(hiv_pos_young, detection=0, prevalence = 0.7)
hiv_pos_aging_ASVs <- core_members(hiv_pos_aging, detection=0, prevalence = 0.7)

hiv_pos_lab_ASVs <- core_members(hiv_pos_lab, detection=0, prevalence = 0.7)
hiv_pos_nonlab_ASVs <- core_members(hiv_pos_nonlab, detection=0, prevalence = 0.7)

hiv_neg_1edu_ASVs <- core_members(hiv_neg_1edu, detection=0, prevalence = 0.7)
hiv_neg_2edu_ASVs <- core_members(hiv_neg_2edu, detection=0, prevalence = 0.7)
hiv_neg_3edu_ASVs <- core_members(hiv_neg_3edu, detection=0, prevalence = 0.7)

hiv_neg_young_ASVs <- core_members(hiv_neg_young, detection=0, prevalence = 0.7)
hiv_neg_aging_ASVs <- core_members(hiv_neg_aging, detection=0, prevalence = 0.7)

hiv_neg_lab_ASVs <- core_members(hiv_neg_lab, detection=0, prevalence = 0.7)
hiv_neg_nonlab_ASVs <- core_members(hiv_neg_nonlab, detection=0, prevalence = 0.7)


# What are these ASVs?

tax_table(prune_taxa(hiv_pos_1edu_ASVs,hiv_final))
tax_table(prune_taxa(hiv_pos_2edu_ASVs,hiv_final))
tax_table(prune_taxa(hiv_pos_3edu_ASVs,hiv_final))
tax_table(prune_taxa(hiv_pos_young_ASVs,hiv_final))
tax_table(prune_taxa(hiv_pos_aging_ASVs,hiv_final))
tax_table(prune_taxa(hiv_pos_lab_ASVs,hiv_final))
tax_table(prune_taxa(hiv_pos_nonlab_ASVs,hiv_final))

tax_table(prune_taxa(hiv_neg_1edu_ASVs,hiv_final))
tax_table(prune_taxa(hiv_neg_2edu_ASVs,hiv_final))
tax_table(prune_taxa(hiv_neg_3edu_ASVs,hiv_final))
tax_table(prune_taxa(hiv_neg_young_ASVs,hiv_final))
tax_table(prune_taxa(hiv_neg_aging_ASVs,hiv_final))
tax_table(prune_taxa(hiv_neg_lab_ASVs,hiv_final))
tax_table(prune_taxa(hiv_neg_nonlab_ASVs,hiv_final))

##HIV POSITIVE##

# can plot those ASVs' relative abundance 
##this is looking funny with everything grouping at tertiary
prune_taxa(hiv_pos_1edu_ASVs,hiv_pos_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Education_Level`, scales ="free")

prune_taxa(hiv_pos_2edu_ASVs,hiv_pos_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Education_Level`, scales ="free")

prune_taxa(hiv_pos_3edu_ASVs,hiv_pos_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Education_Level`, scales ="free")

prune_taxa(hiv_pos_young_ASVs,hiv_pos_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Visit_Age`, scales ="free")

prune_taxa(hiv_pos_aging_ASVs,hiv_pos_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Visit_Age`, scales ="free")

prune_taxa(hiv_pos_lab_ASVs,hiv_pos_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Labor_category`, scales ="free")

prune_taxa(hiv_pos_nonlab_ASVs,hiv_pos_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Labor_category`, scales ="free")


# Notice that in this dataset, there are very few CORE microbiome members. This is common
### What if we wanted to make a Venn diagram of all the ASVs that showed up in each treatment?
hiv_pos_lab_list <- core_members(hiv_pos_lab, detection=0.001, prevalence = 0.50)
hiv_pos_nonlab_list <- core_members(hiv_pos_nonlab, detection=0.001, prevalence = 0.50)

hiv_neg_lab_list <- core_members(hiv_neg_lab, detection=0.001, prevalence = 0.50)
hiv_neg_nonlab_list <- core_members(hiv_neg_nonlab, detection=0.001, prevalence = 0.50)

hiv_pos_aging_list <- core_members(hiv_pos_aging, detection=0.001, prevalence = 0.50)
hiv_pos_young_list <- core_members(hiv_pos_young, detection=0.001, prevalence = 0.50)

hiv_pos_1edu_list <- core_members(hiv_pos_1edu, detection=0.001, prevalence = 0.50)
hiv_pos_2edu_list <- core_members(hiv_pos_2edu, detection=0.001, prevalence = 0.50)
hiv_pos_3edu_list <- core_members(hiv_pos_3edu, detection=0.001, prevalence = 0.50)



hiv_negpos_lab_list_full <- list(HIV_Positive = hiv_pos_lab_list, HIV_Negative = hiv_neg_lab_list)
hiv_lab_nolab_list_full <- list(Labour = hiv_pos_lab_list, No_Labour = hiv_pos_nonlab_list)
hiv_aging_young_list_full <- list(Older_than_45 = hiv_pos_aging_list, Younger_than_45 = hiv_pos_young_list)
hiv_education_full <- list(Primary = hiv_pos_1edu_list, Secondary = hiv_pos_2edu_list, Tertiary = hiv_pos_3edu_list)


library("sf")

# other colours in pallate: "#d7191c", "#fdae61", "#ffffbf", 

custom_colors <-  c("#abd9e9", "#2c7bb6")

hiv_pos_labour_nonlabour_venn <- ggVennDiagram(hiv_lab_nolab_list_full, label_alpha = 0, label_size = 20, set_size = 13) + 
  scale_fill_gradientn(colors = custom_colors) +
  guides(fill = "none")

hiv_pos_young_aging_venn <- ggVennDiagram(hiv_aging_young_list_full, label_alpha = 0, label_size = 20, set_size = 5) + 
  scale_fill_gradientn(colors = custom_colors) +
  guides(fill = "none") 

hiv_pos_education_venn <- ggVennDiagram(hiv_education_full, label_alpha = 0, label_size = 20, set_size = 13) + 
  scale_fill_gradientn(colors = custom_colors) +
  guides(fill = "none") 

ggsave("hiv_pos_labour_nonlabour_venn.png", hiv_pos_labour_nonlabour_venn)
ggsave("hiv_pos_young_aging_venn.png", hiv_pos_young_aging_venn)
ggsave("hiv_pos_education_venn.png", hiv_pos_education_venn)

#################################################################################################
##ISA##

library(tidyverse)
library(phyloseq)
library(indicspecies)

# glom to Genus HIV Positive

hiv_pos_genus <- tax_glom(hiv_pos, "Genus", NArm = FALSE)
hiv_pos_genus_RA <- transform_sample_counts(hiv_pos_genus, fun=function(x) x/sum(x))

#ISA

###HIV POSITIVE###
#Education Level
isa_hiv_pos_edu <- multipatt(t(otu_table(hiv_pos_genus_RA)), cluster = sample_data(hiv_pos_genus_RA)$`Education_Level`)
summary(isa_hiv_pos_edu)
taxtable <- tax_table(hiv_final) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# consider that your table is only going to be resolved up to the genus level, be wary of 
# anything beyond the glomed taxa level
#stat closer to 1, better indicator
isa_hiv_pos_edu_table <- isa_hiv_pos_edu$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) 

#Labour

##Find the file that'll show you which one
isa_hiv_pos_labour <- multipatt(t(otu_table(hiv_pos_genus_RA)), cluster = sample_data(hiv_pos_genus_RA)$`Labor_category`)
summary(isa_hiv_pos_labour)
taxtable <- tax_table(hiv_final) %>% as.data.frame() %>% rownames_to_column(var="ASV")

isa_hiv_pos_labour_table <- isa_hiv_pos_labour$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) 

#Age

isa_hiv_pos_age <- multipatt(t(otu_table(hiv_pos_genus_RA)), cluster = sample_data(hiv_pos_genus_RA)$`Visit_Age`)
summary(isa_hiv_pos_age)
taxtable <- tax_table(hiv_final) %>% as.data.frame() %>% rownames_to_column(var="ASV")

isa_hiv_pos_age_table <- isa_hiv_pos_age$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) 

###HIV NEGATIVE###
#Education Level

#################################################################################################
## DESEQ ##
library(DESeq2)
###HIV POSITIVE###
##AGE##

hiv_pos_plus1 <- transform_sample_counts(hiv_pos, function(x) x+1)
hiv_pos_age_deseq <- phyloseq_to_deseq2(hiv_pos_plus1, ~`Visit_Age`)
DESEQ_hiv_pos_age <- DESeq(hiv_pos_age_deseq)
res <- results(DESEQ_hiv_pos_age, tidy=TRUE, 
               #this will ensure that No is your reference group
               contrast = c("Visit_Age","Older than 45", "Younger than 45"))
View(res)

# Look at results 
custom_colors <- c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6")

## Volcano plot: effect size VS significance
ggplot(res) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
age_vol_plot <- res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

age_vol_plot

ggsave(filename="age_vol_plot.png",age_vol_plot)
# To get table of results
sigASVs <- res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

# Get only asv names
sigASVs_vec <- sigASVs %>%
  pull(ASV)

# Prune phyloseq file
hiv_age_DESeq <- prune_taxa(sigASVs_vec,hiv_pos)
sigASVs <- tax_table(hiv_age_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

sigASVs <- tax_table(hiv_age_DESeq) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs, by = "ASV") %>%
  arrange(log2FoldChange) %>%  # Ensure proper ordering before converting to a factor
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = Genus[order(log2FoldChange)])) %>%  # Explicit ordering
  mutate(Change = factor(ifelse(log2FoldChange > 0, "Positive", "Negative"), 
                         levels = c("Positive", "Negative")))  # Create Change category

# Define colors for positive and negative values
custom_colors <- c("Positive" = "#2c7bb6", "Negative" = "#d7191c")


# Create the bar plot with color mapping
age_bar_plot <- ggplot(sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill=Change), stat="identity") +  
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange - lfcSE, ymax=log2FoldChange + lfcSE)) +
  scale_fill_manual(values=custom_colors) +  
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

# Create the bar plot 
age_bar_plot <- ggplot(sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity") +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange - lfcSE, ymax=log2FoldChange + lfcSE)) +
  scale_fill_manual(values=custom_colors) +  # Assign custom colors
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) 

ggsave("age_bar_plot.png", plot = age_bar_plot)
?ggsave
##LABOR##

hiv_pos_labor_deseq <- phyloseq_to_deseq2(hiv_pos_plus1, ~`Labor_category`)
DESEQ_hiv_pos_labor <- DESeq(hiv_pos_labor_deseq)
res_labor <- results(DESEQ_hiv_pos_labor, tidy=TRUE, 
               #this will ensure that No is your reference group
               contrast = c("Labor_category","Labor", "Non-Labor"))
View(res_labor)

# Look at results 

## Volcano plot: effect size VS significance
ggplot(res_labor) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
labor_vol_plot <- res_labor %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

labor_vol_plot

ggsave(filename="labor_vol_plot.png",labor_vol_plot)
# To get table of results
sigASVs_labor <- res_labor %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_labor)

# Get only asv names
sigASVs_labor_vec <- sigASVs_labor %>%
  pull(ASV)

# Prune phyloseq file
hiv_labor_DESeq <- prune_taxa(sigASVs_labor_vec,hiv_pos)
sigASVs_labor <- tax_table(hiv_labor_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_labor) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

sigASVs_labor <- tax_table(hiv_labor_DESeq) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_labor, by = "ASV") %>%
  arrange(log2FoldChange) %>%  # Ensure proper ordering before converting to a factor
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = Genus[order(log2FoldChange)])) %>%  # Explicit ordering
  mutate(Change = factor(ifelse(log2FoldChange > 0, "Positive", "Negative"), 
                         levels = c("Positive", "Negative")))  # Create Change category

# Define colors for positive and negative values
custom_colors <- c("Positive" = "#2c7bb6", "Negative" = "#d7191c")

# Create the bar plot 
bar_plot_labor <- ggplot(sigASVs_labor) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity") +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange - lfcSE, ymax=log2FoldChange + lfcSE)) +
  scale_fill_manual(values=custom_colors) +  # Assign custom colors
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

bar_plot_labor <- ggplot(sigASVs_labor) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill=Change), stat="identity") +  
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange - lfcSE, ymax=log2FoldChange + lfcSE)) +
  scale_fill_manual(values=custom_colors) +  
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

ggsave("bar_plot_labor.png", plot = bar_plot_labor)

