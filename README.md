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
[Feb 1st](/Meeting_minutes/2025-02-01.md) | [Feb 7th](/Meeting_minutes/2025-02-07.md) | [Feb 11th](/Meeting_minutes/2025-02-11.md) | [Feb 13th](/Meeting_minutes/2025-02-13.md)
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

- Determine ASVs with DADA2
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux_seqs.qza \
  --p-trim-left 0 \
  --p-trunc-len 251 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats.qza








## Finalized Codes


