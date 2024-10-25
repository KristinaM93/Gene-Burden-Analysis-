import pandas
import os
import numpy
from collections import Counter

# Ask inputs
## Paths to input files
tsvPath = "/scratch/RDS-FMH-UNITI-RW/Q1"
samplesPath = "/scratch/RDS-FMH-UNITI-RW/Q1/samples_to_filter.tsv"
genesPath = "all"
## Paths and names to save files
saveDirPath = "/scratch/RDS-FMH-UNITI-RW/Q1"
samples_name = "samples"
gene_name = "genes"


# Read files
## VCF file
### If 'merge' is selected the default TSV will be used, this file contains all the variants for 426 samples
if tsvPath == '/scratch/RDS-FMH-UNITI-RW/Q1':
    print('You have selected the TSV resulting of the VCF merge with all the samples')
    tsvPathUse = '/scratch/RDS-FMH-UNITI-RW/Q1/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.tsv'
    numColsAnnotated = 48 # Number of columns with annotation, not the columns with the genotype, neither columns with frequencies for the original cohort
### User can select the TSV that they desired
else:
    tsvPathUse = tsvPath
    numColsAnnotated = "48"

tsv = pandas.read_csv(tsvPathUse, sep = '\t', index_col = False, dtype = object)
print('The number of positions in the TSV is: ' + str(len(tsv))) # Print information to check if it is correct

# Drop columns with 'cohort', they have the old AC, AF and AN
tsv = tsv[tsv.columns.drop(list(tsv.filter(regex = 'cohort')))]


# Filter by genes
## Keep all the variants (rows), don't filter by gene
if genesPath == 'all':
    print('Any gene will be filtered')
    tsv_filtered_genes = tsv
else:
    ## File with genes
    genes = pandas.read_csv(genesPath, sep = '\t', index_col = False, dtype = object, header = None)
    genes = genes[0].tolist() # Transform to list
    print('The number of genes is: ' + str(len(genes))) # Print information to check if it is correct
    ## Filter
    boolean_series_genes = tsv.SYMBOL.isin(genes)
    tsv_filtered_genes = tsv[boolean_series_genes]


# Filter by samples
## Keep all samples in the TSV
if samplesPath == 'all':
    print('Any sample will be filtered')
    samplesPathUse = 'ALL'
## Keep sporadic Meniere Disease samples, which are saved into a file
#elif samplesPath == 'SMD':
#    print('SMD has been selected')
#    samplesPathUse = '/mnt/d/exomes_vcf_all_files/cohort_20220117/annotated/files2filter/cohort_17112021_SMD_samples.txt'
## Keep familial Meniere Disease samples, which are saved into a file
#elif samplesPath == 'FMD':
#    print('FMD has been selected')
#    samplesPathUse = '/mnt/d/exomes_vcf_all_files/cohort_20220117/annotated/files2filter/cohort_17112021_FMD_samples.txt'
## Keep desired samples
else:
    samplesPathUse = samplesPath
    
## Maintain all the samples
if samplesPathUse == 'ALL':
    tsv_filtered_samples = tsv_filtered_genes
    nSamples = 75
    samples = tsv_filtered_samples.columns[len(tsv_filtered_samples.columns)-nSamples:].tolist()
## Filter by samples
else:
    ## File with samples
    samples = pandas.read_csv(samplesPathUse, sep='\t', index_col=False, dtype=object, header=None)
    samples = samples[0].tolist()  # Transform to list
    print('The number of samples is: ' + str(len(samples)))  # Print information to check if it is correct

    # Order the samples based on their order in the annotated file
    samples_order = [col for col in tsv.columns if col in samples]

    # Filter - keep columns with SNV information and desired samples columns in the correct order
    colsKeep = []
    colsKeep.extend(tsv_filtered_genes.columns.values.tolist()[:numColsAnnotated])
    colsKeep.extend(samples_order)
    tsv_filtered_samples = tsv_filtered_genes[colsKeep]
    nSamples = len(samples_order)

# Convert to float specific columns
## Convert '-' to NaN, '-' cause some errors
tsv_filtered_rows = tsv_filtered_samples.replace('-', numpy.NaN)

colsToFloat = ['CADD_PHRED', 'CADD_RAW', 'gnomADg_AF_nfe', 'gnomADg_AF', 'CSVS_AF']
for col in colsToFloat:
    print(col)
    tsv_filtered_rows[col] = tsv_filtered_rows[col].astype(float)

# Change '|' by '/' in GT(genotype)
for sample in samples:
    tsv_filtered_rows[sample] = tsv_filtered_rows[sample].str.replace('|', '/')

# Allele Count (AC) based on GT
## Create new DF with GT and remove the rest
GT = tsv_filtered_rows[samples]

## Keep 3 first characters
for sample in samples:
    GT[sample] = GT[sample].str[:3]


## Count heterocygous, homogygous, AC (Allele count) and in what samples they are; per row (SNV)
### Parse each genotype for all SNV and count the number of alleles (heterocygous(0/1) = 1 and homocygous(1/1) = 2)
### Sum 1 to number of heterocygous or homocygous and add the sample to the list in each case
### Save which samples are heterocygous or homocygous
#### Convert GT list of lists
GT = GT.replace('1/0', '0/1') # When multiallelics were split, some variants had 1/0 in the genotype, is the same as 0/1

samplesOrder = GT.columns.values.tolist()
GTlist = GT.values.tolist() # Create a list of lists from the DF, to be more efficient
countList = [] # Create empty lists to save the desired data
samples_het = []
samples_hom = []

# For each variant (each list of lists)
for GTl in GTlist:
    # Count each GT: 0/0, 0/1 & 1/1
    count = Counter(GTl)
    countList.append(count)
    # Add to list heterocygous and homocygous samples
    ## Heterocygous
    samples_het_r_index = [i for i, e in enumerate(GTl) if e == '0/1'] # Take the index for heterocygous variants
    samples_het_r = [samplesOrder[i] for i in samples_het_r_index] # Take the name of the sample (column) by the index
    samples_het_r = '|'.join(samples_het_r) # Split samples by '|'
    samples_het.append(samples_het_r)
    ## Homocygous
    samples_hom_r_index = [i for i, e in enumerate(GTl) if e == '1/1'] # Take the index for homocygous variants
    samples_hom_r = [samplesOrder[i] for i in samples_hom_r_index] # Take the name of the sample (column) by the index
    samples_hom_r = '|'.join(samples_hom_r) # Split samples by '|'
    samples_hom.append(samples_hom_r)

### Create DF with count info and change NaN to 0 for the next calculations
countDF = pandas.DataFrame(data = countList)
countDF = countDF.fillna(0.0)

### Add and rename the desired columns in countDF
#### Calculation of AC (allele count), AN (allele number) and AF (allele frequency)
countDF['AC'] = countDF['0/1'] + 2*countDF['1/1'] # Homocygous samples have the variant in the 2 alleles
#countDF['AC'] = 1  # Set AC to a non-zero value
countDF['AN'] = (countDF['0/0'] + countDF['0/1'] + countDF['1/1'])*2 # Two alleles per samples
countDF['AF'] = countDF['AC']/countDF['AN']

countDF = countDF.rename(columns = {'0/1': 'n_het', '1/1': 'n_hom'}) # Number of samples het and hom

countDF['samples_het'] = samples_het
countDF['samples_hom'] = samples_hom

countDF = countDF.drop('0/0', 1) # Not interested in samples without the variant
countDF = countDF.reset_index(drop = True)
print(countDF)


## Add info in countDF to the DF filtered from the input
tsv_filtered = tsv_filtered_rows
tsv_filtered = tsv_filtered.reset_index(drop = True) # Reset index to join correctly
tsv_filtered = tsv_filtered.join(countDF) # Join countDF to tsv_filtered

colsToInt = ['AC', 'AN', 'n_het', 'n_hom'] # Columns to numeric
for col in colsToInt:
    print(col)
    tsv_filtered[col] = tsv_filtered[col].astype(int)


# Remove SNV not present in the subset of patients
tsv_filtered = tsv_filtered.loc[tsv_filtered['AC'] != 0]


# Reorder columns
colsOrder = []
colsOrder.extend(tsv_filtered.columns.values.tolist()[:len(tsv_filtered.columns)-nSamples-len(countDF.columns)])
colsOrder.extend(['AF', 'AC', 'AN', 'n_het', 'n_hom', 'samples_het', 'samples_hom'])
colsOrder.extend(samples)
print(colsOrder)

tsv_filtered = tsv_filtered[colsOrder]
print('The number of positions in the TSV filtered is: ' + str(len(tsv_filtered)))


# SAVE
print('Files will be saved in ' + saveDirPath)
tsvName = os.path.splitext(tsvPathUse.split('/')[-1])[0]

name_save = '.'.join([tsvName, gene_name, samples_name, 'tsv']) # Add the name for the genes and for the samples that the user set as input
name_save_path = os.path.join(saveDirPath, name_save)
tsv_filtered.to_csv(name_save_path, header = True, index = None, sep = '\t', float_format = '%f')
print(name_save + ' saved')

print('PROCESS FINISHED')
