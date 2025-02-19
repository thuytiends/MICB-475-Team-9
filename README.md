# MICB 475: Data Science Research in Microbiology and Immunology

## Timeline
### tentative
<img src="https://github.com/Anmol-Baranwal/Cool-GIFs-For-GitHub/assets/74038190/ff1b5f32-9420-4dde-b2b9-ed2c0aa17459" width="500">
<br><br> 

## Table of Content 🔖
  * [Meeting minutes](#meeting-minutes-and-agenda)
  * [Lab notebook](#lab-notebook)
  * [Finalized codes](#finalized-codes)


## Meeting Minutes and Agenda
### February ⛄
[Feb 1st](/Meeting_minutes/2025-02-01.md) | [Feb 7th](/Meeting_minutes/2025-02-07.md) | [Feb 11th](/Meeting_minutes/2025-02-11.md) | [Feb 13th](/Meeting_minutes/2025-02-13.md) | [Feb 19th](/Meeting_minutes/2025-02-19.md)
### March 🌸

### April 💻

## Lab Notebook 

QIIME2 Analysis

- Create a directory for the project:
mkdir HIV

- Import a dataset using a manifest file:
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path /mnt/datasets/project_2/hiv/hiv_manifest.tsv \
  --output-path demux_seqs.qza

- Create visualization of demultiplexed samples:
qiime demux summarize \
  --i-data demux_seqs.qza \
  --o-visualization demux.qzv

- Running jobs in the background
screen -S Susana_243
screen -S tien-denoise

- Determine ASVs with DADA2
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

- Visualize DADA2 stats
qiime metadata tabulate \
  --m-input-file 243_stats.qza \
  --o-visualization 243_stats.qzv

qiime metadata tabulate \
  --m-input-file stats251.qza \
  --o-visualization stats251.qzv

- Visualize ASVs stats
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


## Finalized Codes


