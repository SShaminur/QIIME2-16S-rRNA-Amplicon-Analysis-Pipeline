#!/bin/bash

# QIIME2 Pipeline Script
# This script performs 16S rRNA amplicon analysis and PICRUSt2 functional prediction

# Exit on error
set -e

# Create output directory structure
mkdir -p q2-picrust2_output ko_metagenome pathabun_exported pathabun_core_metrics_out

echo "Starting QIIME2 Pipeline..."
echo "=========================="

# 1. Data Import and Processing
echo "Step 1: Importing and processing data..."

# 1.1 Import data
echo "1.1 Importing data..."
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest \
  --output-path 1_1_demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

# 1.2 Join paired end reads
echo "1.2 Joining paired end reads..."
qiime vsearch merge-pairs \
  --i-demultiplexed-seqs 1_1_demux.qza \
  --o-merged-sequences 1_2_demux-joined.qza \
  --o-unmerged-sequences 1_2_demux-unmerged.qza

# 1.3 Quality filter
echo "1.3 Quality filtering..."
qiime quality-filter q-score \
  --i-demux 1_2_demux-joined.qza \
  --p-min-quality 20 \
  --o-filtered-sequences 1_3_demux-joined-filtered.qza \
  --o-filter-stats 1_3_demux-joined-filter-stats.qza

# 1.4 Dereplicate sequences
echo "1.4 Dereplicating sequences..."
qiime vsearch dereplicate-sequences \
  --i-sequences 1_2_demux-joined.qza \
  --o-dereplicated-table 1_4_table.qza \
  --o-dereplicated-sequences 1_4_rep-seqs.qza

# 1.5 De novo clustering
echo "1.5 De novo clustering..."
qiime vsearch cluster-features-de-novo \
  --i-table 1_4_table.qza \
  --i-sequences 1_4_rep-seqs.qza \
  --p-perc-identity 0.99 \
  --p-threads 80 \
  --o-clustered-table 1_5_table-dn-99.qza \
  --o-clustered-sequences 1_5_rep-seqs-dn-99.qza

# 1.6 De novo chimera checking
echo "1.6 Chimera checking..."
qiime vsearch uchime-denovo \
  --i-table 1_5_table-dn-99.qza \
  --i-sequences 1_5_rep-seqs-dn-99.qza \
  --output-dir 1_6_uchime-dn-out

# 1.7 Exclude chimeras
echo "1.7 Excluding chimeras..."

# 1.7a Filter table
qiime feature-table filter-features \
  --i-table 1_5_table-dn-99.qza \
  --m-metadata-file 1_6_uchime-dn-out/nonchimeras.qza \
  --o-filtered-table 1_7a_table-dn-99.qza

# 1.7b Filter sequences
qiime feature-table filter-seqs \
  --i-data 1_5_rep-seqs-dn-99.qza \
  --m-metadata-file 1_6_uchime-dn-out/nonchimeras.qza \
  --o-filtered-data 1_7b_rep-seqs-dn-99.qza

# 1.7c Summarize table
qiime feature-table summarize \
  --i-table 1_7a_table-dn-99.qza \
  --o-visualization 1_7a_table-dn-99.qzv

echo "Step 1 completed!"
echo "================="

# 3. Taxonomic Assignment (skip if 3_3_ref-seqs.qza exists)
if [ ! -f "3_3_ref-seqs.qza" ]; then
    echo "Step 3: Taxonomic assignment..."
    echo "3.1-3.2: Download and prepare reference datasets manually"
    echo "Note: Download gg_13_5/gg_13_8_otus.tar.gz from:"
    echo "ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz"
    echo "Extract 99_otus.fasta and 99_otu_taxonomy.txt"
    echo ""
    echo "Please ensure 99_otus.fasta and 99_otu_taxonomy.txt are in the current directory"
    read -p "Press Enter to continue if files are ready, or Ctrl+C to exit..."
    
    # 3.2 Import reference datasets
    echo "3.2 Importing reference datasets..."
    
    # 3.2a Import sequences
    qiime tools import \
      --type 'FeatureData[Sequence]' \
      --input-path 99_otus.fasta \
      --output-path 3_2a_ref-99_otus.qza
    
    # 3.2b Import taxonomy
    qiime tools import \
      --type 'FeatureData[Taxonomy]' \
      --input-format HeaderlessTSVTaxonomyFormat \
      --input-path 99_otu_taxonomy.txt \
      --output-path 3_2b_ref-99_taxonomy.qza
    
    # 3.3 Extract reference reads
    echo "3.3 Extracting reference reads (time-consuming)..."
    qiime feature-classifier extract-reads \
      --i-sequences 3_2a_ref-99_otus.qza \
      --p-f-primer CCTACGGGNGGCWGCAG \
      --p-r-primer GACTACHVGGGTATCTAATCC \
      --p-min-length 300 \
      --p-max-length 500 \
      --o-reads 3_3_ref-seqs.qza
else
    echo "Step 3: 3_3_ref-seqs.qza exists, skipping taxonomic reference preparation..."
fi

# Continue with taxonomic analysis
echo "Continuing with taxonomic analysis..."

# 3.4 Train classifier
echo "3.4 Training classifier..."
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads 3_3_ref-seqs.qza \
  --i-reference-taxonomy 3_2b_ref-99_taxonomy.qza \
  --o-classifier 3_4_classifier.qza

# 3.5 Taxonomic analysis
echo "3.5 Taxonomic classification..."

# 3.5a Classification
qiime feature-classifier classify-sklearn \
  --i-classifier 3_4_classifier.qza \
  --i-reads 1_7b_rep-seqs-dn-99.qza \
  --o-classification 3_5a_taxonomy.qza

# 3.5b Summary
qiime metadata tabulate \
  --m-input-file 3_5a_taxonomy.qza \
  --o-visualization 3_5b_taxonomy.qzv

# 3.6 Bar plot
echo "3.6 Creating taxa bar plots..."
qiime taxa barplot \
  --i-table 1_7a_table-dn-99.qza \
  --i-taxonomy 3_5a_taxonomy.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization 3_6_taxa-bar-plots.qzv

echo "Step 3 completed!"
echo "================="

# 2. PICRUSt2 Preparation
echo "Step 2: Preparing for PICRUSt2..."

# Filter features
echo "Filtering features for PICRUSt2..."
qiime feature-table filter-features \
  --i-table 1_7a_table-dn-99.qza \
  --p-min-samples 1 \
  --o-filtered-table 7_1filteredtable.qza

# Filter sequences
qiime feature-table filter-seqs \
  --i-data 1_7b_rep-seqs-dn-99.qza \
  --i-table 7_1filteredtable.qza \
  --o-filtered-data 7_2filltered_seq.qza

# Filter taxonomy table
qiime taxa filter-table \
  --i-table 7_1filteredtable.qza \
  --i-taxonomy 3_5a_taxonomy.qza \
  --p-include p__ \
  --p-exclude 'p__;,Eukaryota,Chloroplast,Mitochondria,Unassigned,Unclassified' \
  --o-filtered-table 8_1filtered_table.qza

# Filter taxonomy sequences
qiime taxa filter-seqs \
  --i-sequences 7_2filltered_seq.qza \
  --i-taxonomy 3_5a_taxonomy.qza \
  --p-include p__ \
  --p-exclude 'p__;,Eukaryota,Chloroplast,Mitochondria,Unassigned,Unclassified' \
  --o-filtered-sequences 8_2filltered_seq.qza

echo "Step 2 completed!"
echo "================="

# 6. PICRUSt2 Pipeline
echo "Step 6: Running PICRUSt2 pipeline..."
qiime picrust2 full-pipeline \
   --i-table 8_1filtered_table.qza \
   --i-seq 8_2filltered_seq.qza \
   --output-dir q2-picrust2_output \
   --p-placement-tool sepp \
   --p-threads 80 \
   --p-hsp-method mp \
   --p-min-align 0.8 \
   --p-max-nsti 2 \
   --verbose \
   --p-edge-exponent 0

# 7. Summarize PICRUSt2 results
echo "Step 7: Summarizing PICRUSt2 results..."

# Pathway abundance summary
qiime feature-table summarize \
   --i-table q2-picrust2_output/pathway_abundance.qza \
   --o-visualization q2-picrust2_output/pathway_abundance.qzv

# KO metagenome summary
qiime feature-table summarize \
   --i-table q2-picrust2_output/ko_metagenome.qza \
   --o-visualization q2-picrust2_output/ko_metagenome.qzv

# Export KO metagenome
qiime tools export \
   --input-path q2-picrust2_output/ko_metagenome.qza \
   --output-path ko_metagenome

# Convert to TSV
biom convert \
   -i ko_metagenome/feature-table.biom \
   -o ko_metagenome/ko_metagenome.tsv \
   --to-tsv

echo "Step 7 completed!"
echo "================="

# 8. Diversity analysis
echo "Step 8: Diversity analysis..."
echo "Note: Using sampling depth 226702. Adjust if needed."

qiime diversity core-metrics \
   --i-table q2-picrust2_output/pathway_abundance.qza \
   --p-sampling-depth 226702 \
   --m-metadata-file metadata.tsv \
   --output-dir pathabun_core_metrics_out \
   --p-n-jobs 1

# 9. Export pathway abundance
echo "Step 9: Exporting pathway abundance..."
qiime tools export \
   --input-path q2-picrust2_output/pathway_abundance.qza \
   --output-path pathabun_exported

# 10. Convert to TSV
echo "Step 10: Converting to TSV..."
biom convert \
   -i pathabun_exported/feature-table.biom \
   -o pathabun_exported/feature-table.biom.tsv \
   --to-tsv

echo "========================================"
echo "Pipeline completed successfully!"
echo "========================================"
echo ""
echo "Output files created:"
echo "1. Taxonomic results: 3_5a_taxonomy.qza, 3_6_taxa-bar-plots.qzv"
echo "2. PICRUSt2 results in: q2-picrust2_output/"
echo "3. KO metagenome in: ko_metagenome/"
echo "4. Pathway abundance in: pathabun_exported/"
echo "5. Diversity metrics in: pathabun_core_metrics_out/"
echo ""
echo "To view visualizations:"
echo "qiime tools view 3_6_taxa-bar-plots.qzv"
echo "qiime tools view q2-picrust2_output/pathway_abundance.qzv"
echo "qiime tools view q2-picrust2_output/ko_metagenome.qzv"
