# Gene-Burden-Analysis-

# The vcf_annotated_asDF_pick_allele_gene.py processes and annotates a VCF file containing variant data. Here’s a summary of each section:

1. Read the VCF file: Reads the VCF file in two parts—header and data—skipping the initial 3473 rows after the header.
2. Extract Annotations: Isolates variant annotations from the INFO column, which contain information about each variant after "CSQ=". Splits these annotations by commas and creates rows for each.
3. Reformat the DataFrame: Drops unnecessary columns (QUAL, FILTER, INFO, FORMAT), resets the index, and splits annotation columns by the delimiter "|".
4. Set Annotation Column Names: Retrieves annotation field names from the VCF header, matching them to the column format extracted from the ##INFO=<ID=CSQ line in the header.
5. Reorder and Save Data: Merges split annotations back into the main DataFrame, reorders columns, and saves the output as a .tsv file.
6. Handle Ampersands (&): Identifies columns with "CSVS" in their name, replaces cells containing "&" with NaN values, and saves this modified version in another .tsv file.

# The filter_genes_samples_vcfMerge.py processes a large TSV file containing genomic variant data, applying various filters and calculations based on genes and samples provided by the user. Here’s a breakdown of its main steps:

1. Input Handling: It specifies paths for input files (variant data, samples to filter, and genes) and output directories.
2. File Reading: It reads a TSV file, which contains merged VCF data of genomic variants, and optionally filters specific columns.
3. Gene Filtering: If a gene filter file is provided, the code filters variants by the genes listed; otherwise, it retains all genes.
4. Sample Filtering: It filters variants by sample if a sample list is provided, otherwise retaining all samples.
5. Genotype Analysis: For each sample, it converts genotype formats, counts heterozygous (0/1) and homozygous (1/1) genotypes, and calculates allele counts (AC), allele numbers (AN), and allele frequencies (AF).
6. Result Assembly: It merges the calculated genotype information with the filtered TSV data.
7. Output: Finally, it saves the processed and filtered data to a new TSV file with a user-defined naming convention.
