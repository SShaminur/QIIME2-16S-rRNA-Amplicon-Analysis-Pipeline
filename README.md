# QIIME2 16S rRNA Amplicon Analysis Pipeline with PICRUSt2 Functional Prediction

A comprehensive, automated pipeline for processing 16S rRNA amplicon sequencing data from raw FASTQ files to taxonomic classification and functional prediction using QIIME2 and PICRUSt2.

## 🚀 Features

- **End-to-End Processing**: From raw sequencing data to biological insights
- **Paired-End Read Processing**: Merge, quality filter, and denoise paired-end sequences
- **ASV Generation**: De novo clustering at 99% identity with chimera removal
- **Taxonomic Classification**: Using GreenGenes 13_8 reference database
- **Functional Prediction**: PICRUSt2 for metagenome functional content prediction
- **Diversity Analysis**: Core metrics calculation for pathway abundance
- **Multiple Output Formats**: QIIME2 artifacts, visualizations, and TSV exports

## 📋 Pipeline Steps

1. **Data Import & Processing**
   - Import paired-end FASTQ files
   - Join paired-end reads
   - Quality filtering (Q-score ≥20)
   - Dereplication and de novo clustering

2. **Taxonomic Analysis**
   - Reference-based classification using GreenGenes
   - Taxonomic bar plots visualization
   - Chimera filtering and quality control

3. **Functional Prediction (PICRUSt2)**
   - Filtered ASV table preparation
   - Metagenome functional prediction
   - Pathway abundance and KO term analysis

4. **Diversity & Export**
   - Alpha and beta diversity metrics
   - Export to various formats for downstream analysis

## 🛠️ Requirements

- QIIME2 (2023.5 or later recommended)
- PICRUSt2 QIIME2 plugin
- Sufficient computational resources:
  - 32GB+ RAM recommended
  - 80+ CPU cores for optimal performance
  - Adequate storage for intermediate files

## 📁 Input Files Required

- `manifest` - QIIME2 manifest file for FASTQ imports
- `metadata.tsv` - Sample metadata in QIIME2 format
- (Optional) GreenGenes reference files if not pre-downloaded:
  - `99_otus.fasta`
  - `99_otu_taxonomy.txt`

## 🚀 Quick Start

```bash
# Make script executable
chmod +x qiime2_pipeline.sh

# Run the pipeline
./qiime2_pipeline.sh
```

## 📊 Outputs

- **QIIME2 Artifacts (`.qza`)** - Processed data at each step
- **Visualizations (`.qzv`)** - Interactive plots and summaries
- **TSV Exports** - Tab-separated files for external analysis
- **PICRUSt2 Results** - Functional predictions and pathway abundances

## 🔧 Customization

- Adjust `--p-threads` for your CPU cores
- Modify primer sequences for different target regions
- Change clustering identity threshold as needed
- Adjust sampling depth for diversity analysis

## 📝 Citation

If using this pipeline, please cite:
- QIIME2: Bolyen et al., 2019, Nature Biotechnology
- PICRUSt2: Douglas et al., 2020, Nature Biotechnology
- VSEARCH: Rognes et al., 2016, PeerJ
- GreenGenes: DeSantis et al., 2006, Applied and Environmental Microbiology

## ⚠️ Notes

- The script includes checkpoints to skip already-completed steps
- Reference database preparation is time-consuming (~1-2 hours)
- PICRUSt2 analysis requires significant computational resources
- Ensure metadata file matches sample names in the manifest

## 🤝 Contributing

Issues and pull requests are welcome for improvements and bug fixes.
