# Gene-Burden-Analysis-

# The vcf_annotated_asDF_pick_allele_gene.py processes and annotates a VCF file containing variant data. Here’s a summary of each section:

1. Read the VCF file: Reads the VCF file in two parts—header and data—skipping the initial 3473 rows after the header.
2. Extract Annotations: Isolates variant annotations from the INFO column, which contain information about each variant after "CSQ=". Splits these annotations by commas and creates rows for each.
3. Reformat the DataFrame: Drops unnecessary columns (QUAL, FILTER, INFO, FORMAT), resets the index, and splits annotation columns by the delimiter "|".
4. Set Annotation Column Names: Retrieves annotation field names from the VCF header, matching them to the column format extracted from the ##INFO=<ID=CSQ line in the header.
5. Reorder and Save Data: Merges split annotations back into the main DataFrame, reorders columns, and saves the output as a .tsv file.
6. Handle Ampersands (&): Identifies columns with "CSVS" in their name, replaces cells containing "&" with NaN values, and saves this modified version in another .tsv file.
