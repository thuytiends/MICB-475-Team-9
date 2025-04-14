# 1. Load packages
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggplot2)
library(ape)
library(picante)
library(reshape2)

install.packages("ggpubr")  
library(ggpubr)

# 2. Load in the HIV metadata,  OTU table, taxonomy file, and phylogenetic tree 
metafp <- "hiv_metadata.tsv"
meta <- read_delim(metafp, delim="\t")
otufp <- "hiv-feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)
taxfp <- "hiv-taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "hiv-tree.nwk"
phylotree <- read.tree(phylotreefp)

#adjust OTU table (save everything except first column (OTU ID)) to matrix:
otu_mat <- as.matrix(otu[,-1])

# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`

# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

#adjust metadata (Save everything except sampleid as new data frame):
samp_df <- as.data.frame(meta[,-1])

# Make sampleids the rownames
rownames(samp_df)<- meta$'sample-id'

# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

#adjust taxonomy (Convert taxon strings to a table with separate taxa rank columns):
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix()

# Save everything except feature IDs
tax_mat <- tax_mat[,-1]

# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`

# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

#### Create phyloseq object ####
# Merge all into a phyloseq object
hiv <- phyloseq(OTU, SAMP, TAX, phylotree)

# View components of phyloseq object with the following commands
otu_table(hiv)
sample_data(hiv)
tax_table(hiv)
phy_tree(hiv)
otu_table(hiv)
sample_data(hiv)
tax_table(hiv)
phy_tree(hiv)

#### Filter metadata ####
#duplicate phylosesq object
phyloseq_hiv <- hiv_rare

# Remove non-bacterial sequences, if any
hiv_filt <- subset_taxa(hiv,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
         
# Remove samples with less than 100 reads
hiv_final <- prune_samples(sample_sums(hiv_filt)>100, hiv_filt)

#Include only Week 0 patient information by filtering for Initial Visit (Week 0)
filtered_hiv_week0 <- subset_samples(phyloseq_hiv, Visit_Cat == "2nd Visit")
# change labor categories by grouping certain categories under "Non-Labor" and then updating the metadata in filtered_hiv_week0.
labor_rename <- as.data.frame(sample_data(filtered_hiv_week0))

# Ensure Labor_category is a character (not a factor)
labor_rename$Labor_category <- as.character(labor_rename$Labor_category)

# Fix potential whitespace issues (just in case)
labor_rename$Labor_category <- trimws(labor_rename$Labor_category)

# Apply the category renaming
labor_rename$Labor_category <- ifelse(
  labor_rename$Labor_category %in% c("Non_Labor_Houseworker", "Non_Labor_Other", "Unemployed"),
  "Non-Labor",
  "Labor"  # Assign all other values (including NA) as "Labor"
)
# Convert back to a factor 
labor_rename$Labor_category <- factor(labor_rename$Labor_category, levels = c("Labor", "Non-Labor"))

# Check if the renaming worked
table(labor_rename$Labor_category)

# Assign back to phyloseq object
sample_data(filtered_hiv_week0) <- sample_data(labor_rename)

#### Subset into Labor and Non-Labor groups####
# Subset the phyloseq object by HIV_Status 
hiv_pos <- subset_samples(filtered_hiv_week0, HIV_Status == "Positive")
hiv_neg <- subset_samples(filtered_hiv_week0, HIV_Status == "Negative")

# Now subset by Labor_category within the HIV-positive samples (Labor group)
hiv_pos_labor <- subset_samples(hiv_pos, Labor_category == "Labor")
hiv_pos_nonlabor <- subset_samples(hiv_pos, Labor_category == "Non-Labor")
unique(sample_data(hiv_pos_nonlabor)$Labor_category)

# Subset by Labor_category within the HIV-negative samples (Labor group)
hiv_neg_labor <- subset_samples(hiv_neg, Labor_category == "Labor")
hiv_neg_nonlabor <- subset_samples(hiv_neg, Labor_category == "Non-Labor")

# Check if subsets have samples
print(paste("Labor samples:", nsamples(hiv_labor)))
print(paste("Non-Labor samples:", nsamples(hiv_nonlabor)))

####Subset Age groups####
# Convert Visit_Age to numeric (in case it's stored as a character or factor)
sample_data(filtered_hiv_week0)$Visit_Age <- as.numeric(as.character(sample_data(filtered_hiv_week0)$Visit_Age))

Age_groups <- as.data.frame(sample_data(filtered_hiv_week0))

# Categorize into Aging Adults (46-70) and Young Adults (18-45)
Age_groups$Age_Group <- ifelse(
  Age_groups$Visit_Age >= 46 & Age_groups$Visit_Age <= 70, "Aging Adults",
  ifelse(Age_groups$Visit_Age >= 18 & Age_groups$Visit_Age <= 45, "Young Adults", NA)
)

sample_data(filtered_hiv_week0) <- sample_data(Age_groups)

# Check if the grouping worked
table(Age_groups$Age_Group, useNA = "always")

# Subset for Young Adults (18-45 yrs)
hiv_pos_young_adults <- subset_samples(hiv_pos, Age_Group == "Young Adults")
hiv_neg_young_adults <- subset_samples(hiv_neg, Age_Group == "Young Adults")

# Subset for Aging Adults (46-70 yrs)
hiv_pos_aging_adults <- subset_samples(hiv_pos, Age_Group == "Aging Adults")
hiv_neg_aging_adults <- subset_samples(hiv_neg, Age_Group == "Aging Adults")


#Rarefy samples
# rngseed sets a random number. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
# t transposes the table to use rarecurve function
# cex decreases font size
rarecurve(t(as.data.frame(otu_table(hiv_final))), cex=0.1)

hiv_rare <- rarefy_even_depth(hiv_final, rngseed = 1, sample.size = 12082)

save(hiv_final, file="hiv_final.RData")
save(hiv_rare, file="hiv_rare.RData")

#### Load in RData ####
load("hiv_rare.RData")
load("hiv_final.RData")

####Beta diversity Analysis####

##Calculate Weighted UniFrac Distance##
bc_dm_weighted <- distance(filtered_hiv_week0, method = "unifrac", weighted = TRUE) #DELETE MAYBE

#for HIV Status
bc_dm_weighted_pos <- distance(hiv_pos, method = "unifrac", weighted = TRUE) #DELETE MAYBE
bc_dm_weighted_neg <- distance(hiv_neg, method = "unifrac", weighted = TRUE) #DELETE MAYBE

#for HIV Labor and Nonlabor#
bc_dm_weighted_pos_labor <- distance(hiv_pos_labor, method = "unifrac", weighted = TRUE)
bc_dm_weighted_pos_nonlabor <- distance(hiv_pos_nonlabor, method = "unifrac", weighted = TRUE)


bc_dm_weighted_neg_labor <- distance(hiv_neg_labor, method = "unifrac", weighted = TRUE)
bc_dm_weighted_neg_nonlabor <- distance(hiv_neg_nonlabor, method = "unifrac", weighted = TRUE)

#for HIV age
bc_dm_weighted_pos_young <- distance(hiv_pos_young_adults, method = "unifrac", weighted = TRUE)
bc_dm_weighted_neg_young <- distance(hiv_neg_young_adults, method = "unifrac", weighted = TRUE)

bc_dm_weighted_pos_aging <- distance(hiv_pos_aging_adults, method = "unifrac", weighted = TRUE)
bc_dm_weighted_neg_aging <- distance(hiv_neg_aging_adults, method = "unifrac", weighted = TRUE)




#Perform PCoA 
ordination_weighted_pos_labor <- ordinate(hiv_pos_labor, method = "PCoA", distance = bc_dm_weighted_pos_labor)
ordination_weighted_pos_nonlabor <- ordinate(hiv_pos_nonlabor, method = "PCoA", distance = bc_dm_weighted_pos_nonlabor)

ordination_weighted_neg_labor <- ordinate(hiv_neg_labor, method = "PCoA", distance = bc_dm_weighted_neg_labor)
ordination_weighted_neg_nonlabor <- ordinate(hiv_neg_nonlabor, method = "PCoA", distance = bc_dm_weighted_neg_nonlabor)

ordination_weighted_neg <- ordinate(hiv_neg, method = "PCoA", distance = bc_dm_weighted_neg)
ordination_weighted_pos <- ordinate(hiv_pos, method = "PCoA", distance = bc_dm_weighted_pos)


#HIV Age group
ordination_weighted_pos_young <- ordinate(hiv_pos_young_adults, method = "PCoA", distance = bc_dm_weighted_pos_young)
ordination_weighted_neg_young <- ordinate(hiv_neg_young_adults, method = "PCoA", distance = bc_dm_weighted_neg_young)

ordination_weighted_pos_aging <- ordinate(hiv_pos_aging_adults, method = "PCoA", distance = bc_dm_weighted_pos_aging)
ordination_weighted_neg_aging <- ordinate(hiv_neg_aging_adults, method = "PCoA", distance = bc_dm_weighted_neg_aging)



## Create a PCoA plot using ggplot2##
# Create PCoA plot for HIV-Positive samples
pcoa_weighted_plot_HIV_Pos_nonlabor <- plot_ordination(hiv_pos_nonlabor, ordination_weighted_pos_nonlabor, color = "Labor_category") +
  geom_point(size = 3, alpha = 0.8) +
  ggtitle("PCoA - Weighted UniFrac (HIV-Positive)") +
  theme_minimal() +
  scale_color_manual(values = c("#d7191c", "#fdae61"))

# Display the plot
pcoa_weighted_plot_HIV_Pos_nonlabor

# Save the plot
ggsave(filename = "PCoA_HIV_Positive_nonlabor.png", 
       plot = pcoa_weighted_plot_HIV_Pos_nonlabor, 
       height = 4, width = 6)

pcoa_weighted_plot_HIV_Pos_labor <- plot_ordination(hiv_pos_labor, ordination_weighted_pos_labor, color = "Labor_category") +
  geom_point(size = 3, alpha = 0.8) +
  ggtitle("PCoA - Weighted UniFrac (HIV-Positive)") +
  theme_minimal() +
  scale_color_manual(values = c("#d7191c", "#fdae61"))

# Display the plot
pcoa_weighted_plot_HIV_Pos_labor

# Save the plot
ggsave(filename = "PCoA_HIV_Positive_labor.png", 
       plot = pcoa_weighted_plot_HIV_Pos_labor, 
       height = 4, width = 6)

##HIV Negative labor and non-labor PCoA
pcoa_weighted_plot_HIV_Neg_labor <- plot_ordination(hiv_neg_labor, ordination_weighted_neg_labor, color = "Labor_category") +
  geom_point(size = 3, alpha = 0.8) +
  ggtitle("PCoA - Weighted UniFrac (HIV-Negative)") +
  theme_minimal() +
  scale_color_manual(values = c("#d7191c", "#fdae61"))

# Display the plot
pcoa_weighted_plot_HIV_Neg_labor

# Save the plot
ggsave(filename = "PCoA_HIV_Negative_labor.png", 
       plot = pcoa_weighted_plot_HIV_Neg_labor, 
       height = 4, width = 6)

##Non-labor HIV neg##
pcoa_weighted_plot_HIV_Neg_nonlabor <- plot_ordination(hiv_neg_nonlabor, ordination_weighted_neg_nonlabor, color = "Labor_category") +
  geom_point(size = 3, alpha = 0.8) +
  ggtitle("PCoA - Weighted UniFrac (HIV-Negative)") +
  theme_minimal() +
  scale_color_manual(values = c("#d7191c", "#fdae61"))

# Display the plot
pcoa_weighted_plot_HIV_Neg_nonlabor

# Save the plot
ggsave(filename = "PCoA_HIV_Negative_nonlabor.png", 
       plot = pcoa_weighted_plot_HIV_Neg_nonlabor, 
       height = 4, width = 6)



#HIV Status
pcoa_weighted_plot <- plot_ordination(filtered_hiv_week0, ordination_weighted, color = "HIV_Status") +
  geom_point(size = 3, alpha = 0.8) +
  ggtitle("PCoA - Weighted UniFrac (HIV Status)") +
  theme_minimal() +
  scale_color_manual(values = c("#d7191c", "#fdae61"))

# Display the plot
pcoa_weighted_plot


ggsave(filename = "HIVStatus_plot_PCoA.png"
       , pcoa_weighted_plot
       , height=4, width=6)


#HIV Education level#
pcoa_weighted_plot_Education_neg <- plot_ordination(hiv_neg, ordination_weighted_neg, color = "Education_Level") +
  geom_point(size = 3, alpha = 0.8) +
  ggtitle("PCoA - Weighted UniFrac (Education)") +
  theme_minimal()  +
  scale_color_manual(values = c("#d7191c", "#fdae61", "#2c7bb6"))


pcoa_weighted_plot_Education_neg

ggsave(filename = "Education_neg_plot_PCoA.png"
       , pcoa_weighted_plot_Education_neg
       , height=4, width=6)

pcoa_weighted_plot_Education_pos <- plot_ordination(hiv_pos, ordination_weighted_pos, color = "Education_Level") +
  geom_point(size = 3, alpha = 0.8) +
  ggtitle("PCoA - Weighted UniFrac (Education)") +
  theme_minimal()  +
  scale_color_manual(values = c("#d7191c", "#fdae61", "#2c7bb6"))


pcoa_weighted_plot_Education_pos


ggsave(filename = "Education_pos_plot_PCoA.png"
       , pcoa_weighted_plot_Education_pos
       , height=4, width=6)

#HIV labor#

pcoa_weighted_plot_Labor <- plot_ordination(filtered_hiv_week0, ordination_weighted, color = "Labor_category") +
  geom_point(size = 3, alpha = 0.8) +
  ggtitle("PCoA - Weighted UniFrac (Labor Category)") +
  theme_minimal()   +
  scale_color_manual(values = c("#d7191c", "#fdae61", "#2c7bb6"))


pcoa_weighted_plot_Labor

ggsave(filename = "Labor_plot_PCoA.png"
       , pcoa_weighted_plot_Labor
       , height=4, width=6)


#HIV Age groups#
pcoa_weighted_plot_HIV_pos_young <- plot_ordination(hiv_pos_young_adults, ordination_weighted_pos_young, color = "Visit_Age") +
  geom_point(size = 3, alpha = 0.8) +
  ggtitle("PCoA - Weighted UniFrac (HIV-Positive)") +
  theme_minimal() 

# Display the plot
pcoa_weighted_plot_HIV_pos_young

# Save the plot
ggsave(filename = "PCoA_HIV_pos_young.png", 
       plot = pcoa_weighted_plot_HIV_pos_young, 
       height = 4, width = 6)

#HIV Age groups#
pcoa_weighted_plot_HIV_neg_young <- plot_ordination(hiv_neg_young_adults, ordination_weighted_neg_young, color = "Visit_Age") +
  geom_point(size = 3, alpha = 0.8) +
  ggtitle("PCoA - Weighted UniFrac (HIV-Negative)") +
  theme_minimal() 

# Display the plot
pcoa_weighted_plot_HIV_neg_young

# Save the plot
ggsave(filename = "PCoA_HIV_neg_young.png", 
       plot = pcoa_weighted_plot_HIV_neg_young, 
       height = 4, width = 6)

#HIV Age groups#
pcoa_weighted_plot_HIV_pos_aging <- plot_ordination(hiv_pos_aging_adults, ordination_weighted_pos_aging, color = "Visit_Age") +
  geom_point(size = 3, alpha = 0.8) +
  ggtitle("PCoA - Weighted UniFrac (HIV-Positive)") +
  theme_minimal() 

# Display the plot
pcoa_weighted_plot_HIV_pos_aging

# Save the plot
ggsave(filename = "PCoA_HIV_pos_aging.png", 
       plot = pcoa_weighted_plot_HIV_pos_aging, 
       height = 4, width = 6)

#HIV Age groups#
pcoa_weighted_plot_HIV_neg_aging <- plot_ordination(hiv_neg_aging_adults, ordination_weighted_neg_aging, color = "Visit_Age") +
  geom_point(size = 3, alpha = 0.8) +
  ggtitle("PCoA - Weighted UniFrac (HIV-Negative)") +
  theme_minimal() 

# Display the plot
pcoa_weighted_plot_HIV_neg_aging

# Save the plot
ggsave(filename = "PCoA_HIV_neg_aging.png", 
       plot = pcoa_weighted_plot_HIV_neg_aging, 
       height = 4, width = 6)

#### Alpha diversity ####

#HIV Status
HIVstatus_gg_richness <- plot_richness(filtered_hiv_week0, x = "HIV_Status", measures = c("Observed", "Shannon")) +
  xlab("HIV Status") +
  geom_boxplot(aes(fill = HIV_Status)) +  # Ensure fill is mapped
  scale_fill_manual(values = c("#d7191c", "#fdae61"))
HIVstatus_gg_richness

ggsave(filename = "HIVStatus_plot_richness.png"
       , HIVstatus_gg_richness
       , height=4, width=6)

estimate_richness(filtered_hiv_week0)

#Labor Status
LaborStatus_gg_richness <- plot_richness(filtered_hiv_week0, x = "HIV_Status", measures = c("Observed", "Shannon")) +
  xlab("Labor Status") +
  geom_violin(aes(fill = Labor_category), trim = FALSE) +  # Use geom_violin for violin plot
  scale_fill_manual(values = c("#d7191c", "#fdae61"))    # Customize fill colors
  
LaborStatus_gg_richness

ggsave(filename = "LaborStatus_plot_richness.png"
       , LaborStatus_gg_richness
       , height=4, width=6)

#Age
Age_gg_richness <- plot_richness(filtered_hiv_week0, x = "Age_Group", measures = c("Observed", "Shannon")) + 
  xlab("Age Group") +
  geom_violin(aes(fill = Age_Group), alpha = 0.5) +  # Use violin plot
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotates x-axis labels 
  scale_fill_manual(values = c("#d7191c", "#fdae61"))

Age_gg_richness

ggsave(filename = "Visit_Age_plot_richness.png"
       , Age_gg_richness
       , height=4, width=6)


####Education level####
Education_gg_richness <- plot_richness(filtered_hiv_week0, x = "Education_Level", measures = c("Observed", "Shannon")) +
  xlab("Education Level") +
  geom_boxplot(aes(fill = Education_Level)) +  # Ensure fill is mapped
  scale_fill_manual(values = c("#d7191c", "#fdae61", "#2c7bb6")) 

Education_gg_richness

ggsave(filename = "Education_plot_richness.png"
       , Education_gg_richness
       , height=4, width=6)


####calculate Faith's phylogenetic diversity as PD####
phylo_dist <- pd(t(otu_table(filtered_hiv_week0)), phy_tree(filtered_hiv_week0),
                 include.root=F) 
# add PD to metadata table
sample_data(filtered_hiv_week0)$PD <- phylo_dist$PD

#Create HIV-Positive & HIV-Negative Subsets
hiv_pos <- subset_samples(filtered_hiv_week0, HIV_Status == "Positive")
hiv_neg <- subset_samples(filtered_hiv_week0, HIV_Status == "Negative")

# calculate Faith's phylogenetic diversity as PD
phylo_dist_hiv_pos <- pd(t(otu_table(hiv_pos)), phy_tree(hiv_pos),
                 include.root=F) 

phylo_dist_hiv_neg <- pd(t(otu_table(hiv_neg)), phy_tree(hiv_neg),
                 include.root=F) 
# add PD to metadata table
sample_data(hiv_pos)$PD <- phylo_dist$PD
sample_data(hiv_neg)$PD <- phylo_dist$PD


####create subgroups####
#HIV Education Groups
hiv_primary <- subset_samples(sample_data(filtered_hiv_week0), Education_Level == "Primary")
hiv_secondary <- subset_samples(sample_data(filtered_hiv_week0), Education_Level == "Secondary")
hiv_tertiary <- subset_samples(sample_data(filtered_hiv_week0), Education_Level == "Tertiary")

#Age groups
hiv_young <- subset_samples(sample_data(filtered_hiv_week0), Age_Group == "Young Adults")
hiv_aging <- subset_samples(sample_data(filtered_hiv_week0), Age_Group == "Aging Adults")

#Labor
hiv_labor <- subset_samples(sample_data(filtered_hiv_week0), Labor_category == "Labor")
hiv_nonlabor <- subset_samples(sample_data(filtered_hiv_week0), Labor_category == "Non-Labor")


# Ensure proper subsetting while keeping all necessary components
hiv_primary <- prune_samples(sample_data(filtered_hiv_week0)$Education_Level == "Primary", filtered_hiv_week0)
hiv_secondary <- prune_samples(sample_data(filtered_hiv_week0)$Education_Level == "Secondary", filtered_hiv_week0)
hiv_tertiary <- prune_samples(sample_data(filtered_hiv_week0)$Education_Level == "Tertiary", filtered_hiv_week0)

hiv_young <- prune_samples(sample_data(filtered_hiv_week0)$Age_Group == "Young Adults", filtered_hiv_week0)
hiv_aging <- prune_samples(sample_data(filtered_hiv_week0)$Age_Group == "Aging Adults", filtered_hiv_week0)

hiv_labor <- prune_samples(sample_data(filtered_hiv_week0)$Labor_category == "Labor", filtered_hiv_week0)
hiv_nonlabor <- prune_samples(sample_data(filtered_hiv_week0)$Labor_category == "Non-Labor", filtered_hiv_week0)


#Check if Phylogenetic Tree is Present
phy_tree(hiv_primary) 
phy_tree(hiv_secondary)
phy_tree(hiv_tertiary)

phy_tree(hiv_young)
phy_tree(hiv_aging)

phy_tree(hiv_labor)
phy_tree(hiv_nonlabor)

#Assign the Phylogenetic Tree
phy_tree(hiv_primary) <- phy_tree(filtered_hiv_week0)
phy_tree(hiv_secondary) <- phy_tree(filtered_hiv_week0)
phy_tree(hiv_tertiary) <- phy_tree(filtered_hiv_week0)

phy_tree(hiv_young) <- phy_tree(filtered_hiv_week0)
phy_tree(hiv_aging) <- phy_tree(filtered_hiv_week0)

phy_tree(hiv_labor) <- phy_tree(filtered_hiv_week0)
phy_tree(hiv_nonlabor) <- phy_tree(filtered_hiv_week0)

# Faith's PD for each Education Level (including both HIV+ and HIV-)
faith_pd_primary <- pd(t(otu_table(hiv_primary)), phy_tree(hiv_primary), include.root = TRUE)
faith_pd_secondary <- pd(t(otu_table(hiv_secondary)), phy_tree(hiv_secondary), include.root = TRUE)
faith_pd_tertiary <- pd(t(otu_table(hiv_tertiary)), phy_tree(hiv_tertiary), include.root = TRUE)

# Faith's PD for each Age groups (including both HIV+ and HIV-)
faith_pd_young <- pd(t(otu_table(hiv_young)), phy_tree(hiv_young), include.root = TRUE)
faith_pd_aging <- pd(t(otu_table(hiv_aging)), phy_tree(hiv_aging), include.root = TRUE)

# Faith's PD for each Labor groups (including both HIV+ and HIV-)
faith_pd_labor <- pd(t(otu_table(hiv_labor)), phy_tree(hiv_labor), include.root = TRUE)
faith_pd_nonlabor <- pd(t(otu_table(hiv_nonlabor)), phy_tree(hiv_nonlabor), include.root = TRUE)

# Add Faith's PD to the sample metadata for each Education Level group
sample_data(hiv_primary)$PD <- faith_pd_primary$PD
sample_data(hiv_secondary)$PD <- faith_pd_secondary$PD
sample_data(hiv_tertiary)$PD <- faith_pd_tertiary$PD

sample_data(hiv_young)$PD <- faith_pd_young$PD
sample_data(hiv_aging)$PD <- faith_pd_aging$PD

sample_data(hiv_labor)$PD <- faith_pd_labor$PD
sample_data(hiv_nonlabor)$PD <- faith_pd_nonlabor$PD

# Violin plot for Faith's PD 
#Labor groups
#Labor
labor_plot.pd <- ggplot(sample_data(hiv_labor), aes(x = HIV_Status, y = PD, fill = Labor_category)) +
  geom_violin() +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5) +  # Add jittered points
  xlab("HIV Status") +
  ylab("Faith's Phylogenetic Diversity (PD)")  +
  scale_fill_manual(values = c("#fdae61")) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = "center")



# view plot
labor_plot.pd

#save plot
ggsave(filename = "Labor_plot_PD.png"
       , labor_plot.pd
       , height=4, width=6)

#Non-labor
nonlabor_plot.pd <- ggplot(sample_data(hiv_nonlabor), aes(x = HIV_Status, y = PD, fill = Labor_category)) +
  geom_violin() +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5) +  # Add jittered points
  xlab("HIV Status") +
  ylab("Faith's Phylogenetic Diversity (PD)")  +
  scale_fill_manual(values = c("#fdae61")) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = "center")



# view plot
nonlabor_plot.pd

#save plot
ggsave(filename = "Non-Labor_plot_PD.png"
       , nonlabor_plot.pd
       , height=4, width=6)

#Age groups

#Young
young_plot.pd <- ggplot(sample_data(hiv_young), aes(x = HIV_Status, y = PD, fill = Age_Group)) +
  geom_violin() +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5) +  # Add jittered points
  xlab("HIV Status") +
  ylab("Faith's Phylogenetic Diversity (PD)")  +
  scale_fill_manual(values = c("#fdae61")) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = "center")



# view plot
young_plot.pd

#save plot
ggsave(filename = "Young_plot_PD.png"
       , young_plot.pd
       , height=4, width=6)

#Aging
aging_plot.pd <- ggplot(sample_data(hiv_aging), aes(x = HIV_Status, y = PD, fill = Age_Group)) +
  geom_violin() +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5) +  # Add jittered points
  xlab("HIV Status") +
  ylab("Faith's Phylogenetic Diversity (PD)")  +
  scale_fill_manual(values = c("#fdae61")) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = "center")



# view plot
aging_plot.pd

#save plot
ggsave(filename = "Aging_plot_PD.png"
       , aging_plot.pd
       , height=4, width=6)



#Education Level: Primary
primary_plot.pd <- ggplot(sample_data(hiv_primary), aes(x = HIV_Status, y = PD, fill = Education_Level)) +
  geom_violin() +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5) +  # Add jittered points
  xlab("HIV Status") +
  ylab("Faith's Phylogenetic Diversity (PD)")  +
  scale_fill_manual(values = c("#fdae61")) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = "center")

 

# view plot
primary_plot.pd

#save plot
ggsave(filename = "EducationPrimary_plot_PD.png"
       , primary_plot.pd
       , height=4, width=6)

#Education Level: Secondary
secondary_plot.pd <- ggplot(sample_data(hiv_secondary), aes(x = HIV_Status, y = PD, fill = Education_Level)) +
  geom_violin() +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5) +  # Add jittered points
  xlab("HIV Status") +
  ylab("Faith's Phylogenetic Diversity (PD)") +
  scale_fill_manual(values = c("#fdae61")) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = "center")



# view plot
secondary_plot.pd

#save plot
ggsave(filename = "EducationSecondary_plot_PD.png"
       , secondary_plot.pd
       , height=4, width=6)

#Education Level: Tertiary

tertiary_plot.pd <- ggplot(sample_data(hiv_tertiary), aes(x = HIV_Status, y = PD, fill = Education_Level)) +
  geom_violin() +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5) +  # Add jittered points
  xlab("HIV Status") +
  ylab("Faith's Phylogenetic Diversity (PD)") +
  scale_fill_manual(values = c("#fdae61")) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = "center")



# view plot
tertiary_plot.pd

#save plot
ggsave(filename = "EducationTertiary_plot_PD.png"
       , tertiary_plot.pd
       , height=4, width=6)



# Compute Shannon Diversity and Pielou's Evenness
shannon_primary <- estimate_richness(hiv_primary, measures = "Shannon")
shannon_secondary <- estimate_richness(hiv_secondary, measures = "Shannon")
shannon_tertiary <- estimate_richness(hiv_tertiary, measures = "Shannon")

shannon_young <- estimate_richness(hiv_young, measures = "Shannon")
shannon_aging <- estimate_richness(hiv_aging, measures = "Shannon")

shannon_labor <- estimate_richness(hiv_labor, measures = "Shannon")
shannon_nonlabor <- estimate_richness(hiv_nonlabor, measures = "Shannon")


evenness_primary <- shannon_primary$Shannon / log(specnumber(otu_table(hiv_primary)))        
evenness_secondary <- shannon_secondary$Shannon / log(specnumber(otu_table(hiv_secondary)))
evenness_tertiary <- shannon_tertiary$Shannon / log(specnumber(otu_table(hiv_tertiary)))

evenness_young <- shannon_young$Shannon / log(specnumber(otu_table(hiv_young)))
evenness_aging <- shannon_aging$Shannon / log(specnumber(otu_table(hiv_aging)))

evenness_labor <- shannon_labor$Shannon / log(specnumber(otu_table(hiv_labor)))
evenness_nonlabor <- shannon_nonlabor$Shannon / log(specnumber(otu_table(hiv_nonlabor)))

#Extract Sample Data as a Data Frame
shannon_primary_df <- as.data.frame(sample_data(hiv_primary))
shannon_secondary_df <- as.data.frame(sample_data(hiv_secondary))
shannon_tertiary_df <- as.data.frame(sample_data(hiv_tertiary))

evenness_primary_df <- as.data.frame(sample_data(hiv_primary))
evenness_secondary_df <- as.data.frame(sample_data(hiv_secondary))
evenness_tertiary_df <- as.data.frame(sample_data(hiv_tertiary))

evenness_young_df <- as.data.frame(sample_data(hiv_young))
evenness_aging_df <- as.data.frame(sample_data(hiv_aging))

evenness_labor_df <- as.data.frame(sample_data(hiv_labor))
evenness_nonlabor_df <- as.data.frame(sample_data(hiv_nonlabor))

#Include Shannon and Evenness diversity in the metadata
shannon_primary_df$shannon_diversity <- shannon_primary$Shannon
shannon_secondary_df$shannon_diversity <- shannon_secondary$Shannon
shannon_tertiary_df$shannon_diversity <- shannon_tertiary$Shannon

evenness_primary_df$evenness_diversity <- evenness_primary
evenness_secondary_df$evenness_diversity <- evenness_secondary
evenness_tertiary_df$evenness_diversity <- evenness_tertiary

evenness_young_df$evenness_diversity <- evenness_young
evenness_aging_df$evenness_diversity <- evenness_aging

evenness_labor_df$evenness_diversity <- evenness_labor
evenness_nonlabor_df$evenness_diversity <- evenness_nonlabor


str(evenness_primary)
str(evenness_secondary)
str(evenness_tertiary)

####Evenness plots####
# Shannon --> Education Level: Primary
primary_plot.shannon <- ggplot(shannon_primary_df, aes(x = HIV_Status, y = shannon_diversity, fill = Education_Level)) + 
  geom_violin() + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5) +  # Add jittered points
  xlab("HIV Status") + 
  ylab("Shannon Diversity Index") +  
  theme_minimal() +
  scale_fill_manual(values = c("#2c7bb6"))

# view plot
primary_plot.shannon

#save plot
ggsave(filename = "EducationPrimary_plot_Shannon.png"
       , primary_plot.shannon
       , height=4, width=6)


# Shannon -->Education Level: Secondary
secondary_plot.shannon <- ggplot(shannon_secondary_df, aes(x = HIV_Status, y = shannon_diversity, fill = Education_Level)) + 
  geom_violin() + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5) +  # Add jittered points
  xlab("HIV Status") + 
  ylab("Shannon Diversity Index") +  
  theme_minimal() +
  scale_fill_manual(values = c("#2c7bb6"))



# view plot
secondary_plot.shannon

#save plot
ggsave(filename = "EducationSecondary_plot_Shannon.png"
       , primary_plot.shannon
       , height=4, width=6)


# Shannon -->Education Level: Tertiary
tertiary_plot.shannon <- ggplot(shannon_tertiary_df, aes(x = HIV_Status, y = shannon_diversity, fill = Education_Level)) + 
  geom_violin() + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5) +  # Add jittered points
  xlab("HIV Status") + 
  ylab("Shannon Diversity Index") + 
  theme_minimal() +
  scale_fill_manual(values = c("#2c7bb6"))



# view plot
tertiary_plot.shannon

#save plot
ggsave(filename = "EducationTertiary_plot_Shannon.png"
       , primary_plot.shannon
       , height=4, width=6)


# Evenness --> Labor
labor_plot.evenness <- ggplot(evenness_labor_df, aes(x = HIV_Status, y = PD, fill = Labor_category)) +
  geom_violin() +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5) +  # Add jittered points
  xlab("HIV Status") +
  ylab("Pielou's Evenness (J')") +
  theme_minimal() +
  scale_fill_manual(values = c("#ffFfbf")) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = "center")



# view plot
labor_plot.evenness

#save plot
ggsave(filename = "Labor_plot_Evenness.png"
       , labor_plot.evenness
       , height=4, width=6)

# Evenness --> Non-Labor
#Troubleshooting to convert to data frame#
class(evenness_nonlabor_df)
evenness_nonlabor_df <- as.data.frame(sample_data(filtered_hiv_week0))
str(evenness_nonlabor_df)
evenness_nonlabor_df <- as.data.frame(sample_data(filtered_hiv_week0)@.Data)
rownames(evenness_nonlabor_df) <- sample_names(filtered_hiv_week0)  # Keep row names
colnames(evenness_nonlabor_df) <- colnames(sample_data(filtered_hiv_week0))

# Step 1: Perform pairwise Wilcoxon tests
test_results_nonlabor <- compare_means(PD ~ HIV_Status, data = evenness_nonlabor_df, 
                              method = "wilcox.test", paired = FALSE)


# Step 2: Adjust p-values using FDR correction
test_results_nonlabor$p.adj <- p.adjust(test_results_nonlabor$p, method = "fdr")

test_results_nonlabor

# Step 4: Create Violin Plot with Adjusted p-values
nonlabor_plot.evenness <- ggplot(evenness_nonlabor_df, aes(x = HIV_Status, y = PD, fill = Labor_category)) +
  geom_violin() +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5) +  # Add jittered points
  xlab("HIV Status") +
  ylab("Pielou's Evenness (J')") +
  theme_minimal() +
  scale_fill_manual(values = c("#ffFfbf")) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = "center")



# view plot
nonlabor_plot.evenness

#save plot
ggsave(filename = "Non-Labor_plot_Evenness.png"
       , nonlabor_plot.evenness
       , height=4, width=6)


# Evenness --> Young
young_plot.evenness <- ggplot(evenness_young_df, aes(x = HIV_Status, y = PD, fill = Age_Group)) +
  geom_violin() +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5) +  # Add jittered points
  xlab("HIV Status") +
  ylab("Pielou's Evenness (J')") +
  theme_minimal() +
  scale_fill_manual(values = c("#ffFfbf")) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = "center")



# view plot
young_plot.evenness

#save plot
ggsave(filename = "Young_plot_Evenness.png"
       , young_plot.evenness
       , height=4, width=6)

#Aging
aging_plot.evenness <- ggplot(evenness_aging_df, aes(x = HIV_Status, y = PD, fill = Age_Group)) +
  geom_violin() +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5) +  # Add jittered points
  xlab("HIV Status") +
  ylab("Pielou's Evenness (J')") +
  theme_minimal() +
  scale_fill_manual(values = c("#ffFfbf")) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = "center")



# view plot
aging_plot.evenness

#save plot
ggsave(filename = "Aging_plot_Evenness.png"
       , aging_plot.evenness
       , height=4, width=6)



# Evenness --> Education Level: Primary
primary_plot.evenness <- ggplot(evenness_primary_df, aes(x = HIV_Status, y = PD, fill = Education_Level)) +
  geom_violin() +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5) +  # Add jittered points
  xlab("HIV Status") +
  ylab("Pielou's Evenness (J')") +
  theme_minimal() +
  scale_fill_manual(values = c("#ffFfbf")) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = "center")



# view plot
primary_plot.evenness

#save plot
ggsave(filename = "EducationPrimary_plot_Evenness.png"
       , primary_plot.evenness
       , height=4, width=6)

# Evenness --> Education Level: Secondary
secondary_plot.evenness <- ggplot(evenness_secondary_df, aes(x = HIV_Status, y = PD, fill = Education_Level)) +
  geom_violin() +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5) +  # Add jittered points
  xlab("HIV Status") +
  ylab("Pielou's Evenness (J')") +
  theme_minimal() +
  scale_fill_manual(values = c("#ffFfbf")) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = "center")



# view plot
secondary_plot.evenness

#save plot
ggsave(filename = "EducationSecondary_plot_Evenness.png"
       , secondary_plot.evenness
       , height=4, width=6)

# Evenness --> Education Level: Tertiary
tertiary_plot.evenness <- ggplot(evenness_tertiary_df, aes(x = HIV_Status, y = PD, fill = Education_Level)) +
  geom_violin() +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5) +  # Add jittered points
  xlab("HIV Status") +
  ylab("Pielou's Evenness (J')") +
  theme_minimal() +
  scale_fill_manual(values = c("#ffFfbf")) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = "center")



# view plot
tertiary_plot.evenness

#save plot
ggsave(filename = "EducationTertiary_plot_Evenness.png"
       , tertiary_plot.evenness
       , height=4, width=6)


####PERMANOVA####
dm_unifrac <- UniFrac(filtered_hiv_week0, weighted=TRUE) # Weighted UniFrac

#Labor
#Generate pcoa plot Permanova
ordination_weighted_laborcat <- ordinate(filtered_hiv_week0, method = "PCoA", distance = dm_unifrac_labor)
PERMANOVA_laborcat <-plot_ordination(filtered_hiv_week0, ordination_weighted_laborcat, color= "Labor_category") +
  stat_ellipse(type = "norm")

ggsave(filename = "PERMANOVA_plot_laborcat.png"
       , PERMANOVA_laborcat
       , height=4, width=6)
#Convert to dataframe
metadata_df <- as.data.frame(labor_rename)
metadata_df <- as.data.frame(sample_data(filtered_hiv_week0)@.Data)
rownames(metadata_df) <- sample_names(filtered_hiv_week0)  # Keep row names
colnames(metadata_df) <- colnames(sample_data(filtered_hiv_week0))


test_results <- compare_means(Age_Group ~ HIV_Status, data = metadata_df, 
                              method = "wilcox.test", paired = FALSE)

# FDR correction
test_results$p.adj <- p.adjust(test_results$p, method = "fdr")


class(metadata_df) 
#PERMANOVA Analysis
permanova_result <- adonis2(dm_unifrac_labor ~ Labor_category, 
                            data = metadata_df, 
                            permutations = 999)
permanova_result #P=  0.478 --> not significant

#Education
ordination_weighted_education <- ordinate(filtered_hiv_week0, method = "PCoA", distance = dm_unifrac)
#Generate pcoa plot Permanova
PERMANOVA_education <-plot_ordination(filtered_hiv_week0, ordination_weighted_laborcat, color= "Education_Level") +
  stat_ellipse(type = "norm")

PERMANOVA_education

ggsave(filename = "PERMANOVA_plot_education.png"
       , PERMANOVA_education
       , height=4, width=6)

#PERMANOVA Analysis
permanova_result_education <- adonis2(dm_unifrac ~ Education_Level, 
                            data = metadata_df, 
                            permutations = 999)
permanova_result_education #p=0.639 not significant

#HIV Status
permanova_result_HIVStatus <- adonis2(dm_unifrac ~ HIV_Status, 
                                      data = metadata_df, 
                                      permutations = 999)
permanova_result_HIVStatus #p= 0.005

#Generate pcoa plot Permanova
PERMANOVA_HIVStatus <-plot_ordination(filtered_hiv_week0, ordination_weighted_laborcat, color= "HIV_Status") +
  stat_ellipse(type = "norm")

PERMANOVA_HIVStatus

ggsave(filename = "PERMANOVA_plot_HIVStatus.png"
       , PERMANOVA_HIVStatus
       , height=4, width=6)

#Age groups
#Convert to dataframe
Agemetadata_df <- as.data.frame(Age_groups)
Agemetadata_df <- as.data.frame(sample_data(filtered_hiv_week0)@.Data)
rownames(Agemetadata_df) <- sample_names(filtered_hiv_week0)  # Keep row names
colnames(Agemetadata_df) <- colnames(sample_data(filtered_hiv_week0))

class(Agemetadata_df) 

#PERMANOVA Analysis
permanova_result_agegrps <- adonis2(dm_unifrac ~ Age_Group, 
                                      data = Agemetadata_df, 
                                      permutations = 999)
permanova_result_agegrps #p= 0.029

#Generate pcoa plot Permanova
PERMANOVA_AgeGroups <-plot_ordination(filtered_hiv_week0, ordination_weighted_laborcat, color= "Age_Group") +
  stat_ellipse(type = "norm")

PERMANOVA_AgeGroups

ggsave(filename = "PERMANOVA_plot_AgeGroups.png"
       , PERMANOVA_AgeGroups
       , height=4, width=6)



####DESeq####
library(DESeq2)
library(phyloseq)

# Install DESeq2 from Bioconductor
BiocManager::install("DESeq2")

#Age
hiv_deseq_agegrp <- phyloseq_to_deseq2(hiv_pos, ~`Age_Group`)
hiv_agegrp_plus1 <- transform_sample_counts(hiv_pos, function(x) x+1)

hiv_deseq_agegrp <- phyloseq_to_deseq2(hiv_agegrp_plus1, ~`Age_Group`)
DESEQ_hiv_agegrp <- DESeq(hiv_deseq_agegrp)

results_hiv_agegrp <- results(DESEQ_hiv_agegrp, tidy=TRUE, 
                              #this will ensure that No is your reference group
                              contrast = c("Age_Group","Young Adults","Aging Adults"))

View(results_hiv_agegrp)
colnames(colData(DESEQ_hiv_labor))
#Education
hiv_deseq_Education <- phyloseq_to_deseq2(hiv_pos, ~`Education_Level`)
hiv_Education_plus1 <- transform_sample_counts(hiv_pos, function(x) x+1)

hiv_deseq_Education <- phyloseq_to_deseq2(hiv_Education_plus1, ~`Education_Level`)
DESEQ_hiv_Education <- DESeq(hiv_deseq_Education)


results_hiv_Education <- results(DESEQ_hiv_Education, tidy=TRUE, 
                              #this will ensure that No is your reference group
                              contrast = c("Education_Level","Primary","Secondary"))

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
#Age
volcanoplot_hiv_agegrp <- ggplot(results_hiv_agegrp) +   
  geom_point(aes(x = log2FoldChange, y = -log10(padj)), color = "#2c7bb6") +  # Use color instead of fill
  theme(
    axis.text = element_text(size = 16),      # Axis text size
    axis.title = element_text(size = 18)      # Axis title size
  )
   

volcanoplot_hiv_agegrp
ggsave(filename = "DESeq_Volcanoplot_AgeGroups.png"
       , volcanoplot_hiv_agegrp_hivneg
       , height=4, width=6)


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
  theme_bw() +  
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

# Sort by absolute log2FoldChange and get only the top 10 ASVs
top10_sigASVs <- sigASVs %>%
  arrange(desc(abs(log2FoldChange))) %>% 
  slice_head(n = 10)  # Select only the top 10

View(sigASVs)
nrow(sigASVs)
# Get only asv names
sigASVs_vec <- sigASVs %>%
  pull(ASV)


sigASVs_vec <- top10_sigASVs %>% 
  pull(ASV)
# Prune phyloseq file
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
ggsave(filename = "DESeqsig_hiv_Education_sec_ter_hivpos.png", plot = DESeqsig_hiv_Education, width = 12, height = 8, dpi = 300)

####DESeq2 Pathway Analysis####

# Create a list of all the packages you need to install
pkgs <- c("ALDEx2", "SummarizedExperiment", "Biobase", "devtools", 
          "ComplexHeatmap", "BiocGenerics", "metagenomeSeq", 
          "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2")

# Use the above list to install all the packages using a for loop
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

install.packages("/Users/susanacw/Downloads/ggpicrust2_1.7.3.tar.gz", repos = NULL, type = "source")

# Load all necessary libraries
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(DESeq2)
library(ggh4x)

#### Import files and preparing tables ####
#Importing the pathway PICrsut2
abundance_file <- "pathway_abundance.tsv"
abundance_data <- read_delim("pathway_abundance.tsv", delim = "\t", col_names = TRUE, trim_ws = TRUE, skip = 1)
abundance_data  =as.data.frame(abundance_data)
view(abundance_data)


#Example Looking at subject number
#If you have multiple variants, filter your metadata to include only 2 at a time

#Remove NAs for your column of interest in this case subject

metadata = metadata[metadata$Visit_Cat == "2nd Visit",]
metadata = metadata[metadata$HIV_Status == "Positive",]

# Filter the metadata to include only Primary and Secondary Education Levels
metadata_filtered = metadata[metadata$Education_Level %in% c("Primary", "Secondary"),]

metadata$Labor_category <- ifelse(metadata$Labor_category %in% c("Non_Labor_Houseworker","Non_Labor_Other","Unemployed"), "Non-Labor", "Labor")

metadata_filtered = metadata[metadata$Labor_category %in% c("Non-Labor", "Labor"),]

metadata$Visit_Age <- ifelse(metadata$Visit_Age > 45, "Older than 45", "Younger than 45")

metadata = metadata[!is.na(metadata$Visit_Age),]
metadata = metadata[!is.na(metadata$Education_Level),]
metadata_filtered = metadata[!is.na(metadata$Labor_category),]
# Get sample names corresponding to the filtered metadata
sample_names = metadata_filtered$`sample-id`

# Add pathway to the sample names if it's not already present
sample_names = append(sample_names, "pathway")

if (!"pathway" %in% colnames(abundance_data)) {
  colnames(abundance_data)[1] <- "pathway"
} 

# Filtering the abundance table to only include samples that are in the filtered metadata
abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% sample_names]

# Remove individuals with no data that caused a problem for pathways_daa()
abundance_data_filtered = abundance_data_filtered[, colSums(abundance_data_filtered != 0) > 0]

# Ensuring the rownames for the abundance_data_filtered is empty (required for their functions to run)
rownames(abundance_data_filtered) = NULL

# Verify samples in metadata match samples in the filtered abundance data
abun_samples = rownames(t(abundance_data_filtered[,-1]))  # Getting a list of the sample names in the newly filtered abundance data

metadata_filtered = metadata_filtered[metadata_filtered$`sample-id` %in% abun_samples,]  # Making sure the filtered metadata only includes these samples

#Perform pathway DAA using DESEQ2 method
abundance_daa_results_df <- pathway_daa(abundance = abundance_data_filtered %>% column_to_rownames("pathway"), 
                                        metadata = metadata_filtered, group = "Labor_category", daa_method = "DESeq2")

DESeq2_metadata = DESeq2_metadata[!is.na(DESeq2_metadata$Labor_category), ]
colnames(DESeq2_metadata)
unique(metadata_filtered$Visit_Age)
str(metadata_filtered$Labor_category)
colnames(metadata_filtered)


# Annotate MetaCyc pathway so they are more descriptive
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                       daa_results_df = abundance_daa_results_df, ko_to_kegg = FALSE)

# Filter p-values to only significant ones
feature_with_p_0.05 <- abundance_daa_results_df %>% filter(p_values < 0.05)

#Changing the pathway column to description for the results 
feature_desc = inner_join(feature_with_p_0.05,metacyc_daa_annotated_results_df, by = "feature")
feature_desc$feature = feature_desc$description
feature_desc = feature_desc[,c(1:7)]
colnames(feature_desc) = colnames(feature_with_p_0.05)

#Changing the pathway column to description for the abundance table
abundance = abundance_data_filtered %>% filter(pathway %in% feature_with_p_0.05$feature)
colnames(abundance)[1] = "feature"
abundance_desc = inner_join(abundance,metacyc_daa_annotated_results_df, by = "feature")
abundance_desc$feature = abundance_desc$description
#this line will change for each dataset. 34 represents the number of samples in the filtered abundance table
abundance_desc = abundance_desc[,-c(106:ncol(abundance_desc))] 

str(abundance_desc)

view(abundance_desc)
ncol(abundance_data_filtered)
dim(abundance_desc)

#Check for duplicate "feature" names
sum(duplicated(abundance_desc$feature))

#Identify the duplicated features
abundance_desc$feature[duplicated(abundance_desc$feature)]

#Merge duplicate pathways (Summing their values)
abundance_desc <- abundance_desc %>%
  group_by(feature) %>%
  summarise(across(everything(), sum))

str(abundance_desc)

# Generate a heatmap
Education_heatmap_pri_sec <-pathway_heatmap(abundance = abundance_desc %>% column_to_rownames("feature"), metadata = metadata_filtered, group = "Education_Level") 

Education_heatmap_pri_ter <-pathway_heatmap(abundance = abundance_desc %>% column_to_rownames("feature"), metadata = metadata_filtered, group = "Education_Level") 

Education_heatmap_sec_ter <-pathway_heatmap(abundance = abundance_desc %>% column_to_rownames("feature"), metadata = metadata_filtered, group = "Education_Level") 

Labor_heatmap <-pathway_heatmap(abundance = abundance_desc %>% column_to_rownames("feature"), metadata = metadata_filtered, group = "Labor_category") 

Age_heatmap <-pathway_heatmap(abundance = abundance_desc %>% column_to_rownames("feature"), metadata = metadata_filtered, group = "Visit_Age") 

Education_heatmap + 
  scale_fill_gradientn(colors = c("#d7191c", "#fdae61", "#2c7bb6"))

print(Education_heatmap)
ggsave("Education_heatmap_pri_sec.png", plot = Education_heatmap_pri_sec, width = 15, height = 8)

ggsave("Education_heatmap_pri_ter.png", plot = Education_heatmap_pri_ter, width = 15, height = 8)

ggsave("EEducation_heatmap_sec_ter.png", plot = Education_heatmap_sec_ter, width = 15, height = 8)

ggsave("Labor_heatmap.png", plot = Labor_heatmap, width = 15, height = 8)

#top 10 labor
Labor_heatmap_top10 <- pathway_heatmap(abundance = top_features %>% column_to_rownames("feature"), 
                                                   metadata = metadata_filtered, 
                                                   group = "Labor_category")
ggsave("Labor_heatmap_top10.png", plot = Labor_heatmap_top10, width = 20, height = 8)

Labor_heatmap_top10
#top 10 for primary and secondary
top_features <- head(abundance_desc, 10)  # Take top 10 pathways
Education_heatmap_top10_pri_sec <- pathway_heatmap(abundance = top_features %>% column_to_rownames("feature"), 
                                     metadata = metadata_filtered, 
                                     group = "Education_Level")
Education_heatmap_top10_pri_sec
ggsave("Education_heatmap_top10_pri_sec.png", plot = Education_heatmap_top10_pri_sec, width = 20, height = 8)

#top 10 primary and tertiary
Education_heatmap_top10_pri_ter <- pathway_heatmap(abundance = top_features %>% column_to_rownames("feature"), 
                                                   metadata = metadata_filtered, 
                                                   group = "Education_Level")

ggsave("Education_heatmap_top10_pri_ter.png", plot = Education_heatmap_top10_pri_ter, width = 20, height = 8)


#top 10 for secondary and tertiary
Education_heatmap_top10_sec_ter <- pathway_heatmap(abundance = top_features %>% column_to_rownames("feature"), 
                                                   metadata = metadata_filtered, 
                                                   group = "Education_Level")



ggsave("Education_heatmap_top10_sec_ter.png", plot = Education_heatmap_top10_sec_ter, width = 20, height = 8)


# Generate pathway PCA plot
Education_pathway_pca <- pathway_pca(abundance = abundance_data_filtered %>% column_to_rownames("pathway"),
                                     metadata = metadata_filtered, group = "Education_Level")




# Generating a bar plot representing log2FC from the custom deseq2 function

# Go to the Deseq2 function script and update the metadata category of interest

# Lead the function in
source("DESeq2_function.R")  

##Education pathways analysis##
# Run the function on your own data
res =  DEseq2_function(abundance_data_filtered, metadata_filtered, "Labor_category")
res$feature =rownames(res)
res_desc = inner_join(res,metacyc_daa_annotated_results_df, by = "feature")
res_desc = res_desc[, -c(8:13)]

# Filter to only include significant pathways
sig_res = res_desc %>%
  filter(pvalue < 0.05, abs(log2FoldChange) >= 2)  # Added log2FC filter

# Order by log2FoldChange
sig_res <- sig_res[order(sig_res$log2FoldChange),]


# You can also filter by Log2fold change

ggplot(data = sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), 
                           x = log2FoldChange, fill = pvalue)) +
  geom_bar(stat = "identity") +  
  theme_bw() +
  labs(x = "Log2 Fold Change", y = "Pathways")


library(ggrepel)
# Ensure log2FoldChange and pvalue are numeric
res_desc <- res_desc %>%
  mutate(log2FoldChange = as.numeric(log2FoldChange),
         pvalue = as.numeric(pvalue))

# Add significance column for coloring points
res_desc <- res_desc %>%
  mutate(Significance = case_when(
    pvalue < 0.05 & abs(log2FoldChange) >= 2 ~ "Significant",
    pvalue < 0.05 & abs(log2FoldChange) < 2 ~ "Moderate",
    TRUE ~ "Not Significant"
  ))

# Volcano plot for pathways
Education_volcanoplot_pathways <-ggplot(res_desc, aes(x = log2FoldChange, y = -log10(pvalue), color = Significance, label = description)) +
  geom_point(alpha = 0.7, size = 3) +  # Pathway points
  scale_color_manual(values = c("Significant" = "red", "Moderate" = "orange", "Not Significant" = "gray")) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black") +  # Log2FC cutoff
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # P-value cutoff
  geom_text_repel(data = res_desc %>% filter(pvalue < 0.01 & abs(log2FoldChange) > 2.5),
                  aes(label = description), size = 5, max.overlaps = 10, nudge_x = 3, hjust = 0) +  # Labeling most significant pathways
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Pathway Expression",
       x = "Log2 Fold Change (Pathway Enrichment)",
       y = "-Log10(p-value)")

Education_volcanoplot_pathways
ggsave("Education_volcanoplot_pathways.png", plot = Education_volcanoplot_pathways, width = 20, height = 8)

##Labor pathways analysis##
# Run the function on your own data
res =  DEseq2_function(abundance_data_filtered, metadata_filtered, "Labor_category")
res$feature =rownames(res)
res_desc = inner_join(res,metacyc_daa_annotated_results_df, by = "feature")
res_desc = res_desc[, -c(8:13)]

# Filter to only include significant pathways
sig_res = res_desc %>%
  filter(pvalue < 0.05, abs(log2FoldChange) >= 2)  # Added log2FC filter

# Order by log2FoldChange
sig_res <- sig_res[order(sig_res$log2FoldChange),]


# You can also filter by Log2fold change

ggplot(data = sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), 
                           x = log2FoldChange, fill = pvalue)) +
  geom_bar(stat = "identity") +  
  theme_bw() +
  labs(x = "Log2 Fold Change", y = "Pathways")


library(ggrepel)
# Ensure log2FoldChange and pvalue are numeric
res_desc <- res_desc %>%
  mutate(log2FoldChange = as.numeric(log2FoldChange),
         pvalue = as.numeric(pvalue))

# Add significance column for coloring points
res_desc <- res_desc %>%
  mutate(Significance = case_when(
    pvalue < 0.05 & abs(log2FoldChange) >= 2 ~ "Significant",
    pvalue < 0.05 & abs(log2FoldChange) < 2 ~ "Moderate",
    TRUE ~ "Not Significant"
  ))

# Volcano plot for pathways
Education_volcanoplot_pathways <-ggplot(res_desc, aes(x = log2FoldChange, y = -log10(pvalue), color = Significance, label = description)) +
  geom_point(alpha = 0.7, size = 3) +  # Pathway points
  scale_color_manual(values = c("Significant" = "red", "Moderate" = "orange", "Not Significant" = "gray")) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black") +  # Log2FC cutoff
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # P-value cutoff
  geom_text_repel(data = res_desc %>% filter(pvalue < 0.01 & abs(log2FoldChange) > 2.5),
                  aes(label = description), size = 5, max.overlaps = 10, nudge_x = 3, hjust = 0) +  # Labeling most significant pathways
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Pathway Expression",
       x = "Log2 Fold Change (Pathway Enrichment)",
       y = "-Log10(p-value)")

Education_volcanoplot_pathways
ggsave("Education_volcanoplot_pathways.png", plot = Education_volcanoplot_pathways, width = 20, height = 8)
