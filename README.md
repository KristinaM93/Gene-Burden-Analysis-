# Gene Burden Analysis
This repository provides scripts for gene burden analysis on genetic variant data, including annotation, filtering, and merging tasks.

# Contents
1. vcf_annotated_asDF_pick_allele_gene.py: Processes annotated VCFs to extract specific allele and gene information as a DataFrame for further analysis.
2. filter_genes_samples_vcfMerge.py: Merges gene and sample data from VCF files.
3. GBA.AB_GQ_DP.filterCSVS.py: Filters CSV files based on specific metrics for quality control.
4. PBS Scripts (vcf.pbs, filter.pbs, gba.pbs): Scripts for submitting jobs to a PBS cluster.
5. Data files (samples_to_filter.tsv, type_analyses_AF005.tsv): Example data for filtering and analysis.
   
# Usage
1. Modify .py scripts as needed for your data.
2. Submit jobs using the PBS scripts on a cluster environment.
