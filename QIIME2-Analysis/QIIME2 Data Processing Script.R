QIIME2 Analysis

# Create a directory for the project:
mkdir HIV

# Import a dataset using a manifest file:
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path /mnt/datasets/project_2/hiv/hiv_manifest.tsv \
  --output-path demux_seqs.qza

# Create visualization of demultiplexed samples:
qiime demux summarize \
  --i-data demux_seqs.qza \
  --o-visualization demux.qzv

# Running jobs in the background
screen -S Susana_243
screen -S tien-denoise

# Determine ASVs with DADA2
Truncating with 251 nts:

qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux_seqs.qza \
  --p-trim-left 0 \
  --p-trunc-len 251 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats.qza

Truncating with 243 nts:
nohup qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux_seqs.qza \
  --p-trim-left 0 \
  --p-trunc-len 243 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats.qza > dada2_output.log 2>&1 &

# Visualize DADA2 stats
qiime metadata tabulate \
  --m-input-file 243_stats.qza \
  --o-visualization 243_stats.qzv

qiime metadata tabulate \
  --m-input-file stats251.qza \
  --o-visualization stats251.qzv

# Convert files to visualize ASVs stats:
qiime feature-table summarize \
  --i-table 243_table.qza \
  --o-visualization 243_table.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/hiv/hiv_metadata.tsv

  qiime feature-table summarize \
  --i-table table251.qza \
  --o-visualization table251.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/hiv/hiv_metadata.tsv


qiime feature-table tabulate-seqs \
  --i-data 243_rep-seqs.qza \
  --o-visualization 243_rep-seqs.qzv

  qiime feature-table tabulate-seqs \
  --i-data rep-seqs251.qza \
  --o-visualization rep-seqs251.qzv


# Train classifier:

qiime feature-classifier extract-reads \ 
  --i-sequences /mnt/datasets/silva_ref_files/silva-138-99-seqs.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-trunc-len 243 \
  --o-reads ref-seqs-trimmed.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs-trimmed.qza \
  --i-reference-taxonomy /mnt/datasets/silva_ref_files/silva-138-99-tax.qza \
  --o-classifier trained_243nt_silva_classifier.qza

# Assign taxonomy to your reads (rep-seqs.qza) with the classifier generated:

qiime feature-classifier classify-sklearn \
  --i-classifier trained_243nt_silva_classifier.qza \
  --i-reads 243_rep-seqs.qza \
  --o-classification 243_taxonomy.qza

qiime metadata tabulate \
  --m-input-file 243_taxonomy.qza \
  --o-visualization 243_taxonomy.qzv

qiime taxa barplot \
  --i-table 243_table.qza \
  --i-taxonomy 243_taxonomy.qza \
  --m-metadata-file /mnt/datasets/project_2/hiv/hiv_metadata.tsv \
  --o-visualization taxa-bar-plots.qzv

# Taxonomy-based filtering
qiime taxa filter-table \
  --i-table 243_table.qza \
  --i-taxonomy 243_taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table table-no-mitochondria-no-chloroplast.qza

qiime feature-table summarize \
  --i-table table-no-mitochondria-no-chloroplast.qza \
  --o-visualization CORRECT_table-no-mitochondria-no-chloroplast.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/hiv/hiv_metadata.tsv

# Generate a tree for phylogenetic diversity analyses:

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences 243_rep-seqs.qza \
  --o-alignment 243_aligned-rep-seqs.qza \
  --o-masked-alignment 243_masked-aligned-rep-seqs.qza \
  --o-tree 243_unrooted-tree.qza \
  --o-rooted-tree 243_rooted-tree.qza 

# Alpha-rarefaction curve:

qiime diversity alpha-rarefaction \
  --i-table 243_table.qza \
  --i-phylogeny 243_rooted-tree.qza \
  --p-max-depth 89000 \
  --m-metadata-file /mnt/datasets/project_2/hiv/hiv_metadata.tsv \
  --o-visualization alpha-rarefaction.qzv


# Determining alpha- and beta-diversity metrics:

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny 243_rooted-tree.qza \
  --i-table table-no-mitochondria-no-chloroplast.qza \
  --p-sampling-depth 12082 \
  --m-metadata-file /mnt/datasets/project_2/hiv/hiv_metadata.tsv \
  --output-dir 12082-core-metrics-results

# Calculate alpha-group-significance

To visualize faith-pd-group-significance:

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file /mnt/datasets/project_2/hiv/hiv_metadata.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

To visualize evenness-group-significance:

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file /mnt/datasets/project_2/hiv/hiv_metadata.tsv \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

#Export to RStudio

#export taxonomy
qiime tools export \
--input-path /HIV/243_taxonomy.qza \
--output-path taxonomy_export 

#export feature-table.txt
qiime tools export \
--input-path /HIV/243_table.qza \
--output-path table_export 

#export rooted tree 
qiime tools export \
--input-path /HIV/243_rooted-tree.qza \
--output-path rooted_tree_expor

####PICRUSt2 analysis####

#Filter the table.qza to remove all the features with 5 or lower counts
qiime feature-table filter-features \
--i-table 243_table.qza \
--p-min-frequency 5 \
--o-filtered-table 243_feature-frequency-filtered-table.qza

qiime tools export \
  --input-path 243_feature-frequency-filtered-table.qza \
  --output-path 243_filtered_exported_table


biom convert \
  -i 243_filtered_exported_table/feature-table.biom \
  -o 243_filtered_exported_table/feature_table.tsv \
  --to-tsv


qiime picrust2 full-pipeline \
  --i-table 243_feature-frequency-filtered-table.qza \
  --i-seq 243_rep-seqs.qza \
  --output-dir q2-picrust2_output \
  --p-placement-tool sepp \
  --p-hsp-method pic \
  --p-max-nsti 2 \
  --verbose

# Convert the output files to human readable files

qiime tools export \
  --input-path HIV/q2-picrust2_output/pathway_abundance.qza \
  --output-path pathabun_exported

biom convert \
   -i pathabun_exported/feature-table.biom \
   -o pathabun_exported/pathway_abundance.tsv \
   --to-tsv



