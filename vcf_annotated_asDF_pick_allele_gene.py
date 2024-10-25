import pandas
import os
import numpy
from collections import Counter

dirPath = '/scratch/RDS-FMH-UNITI-RW/Q1_UNITI'
vcfPath = os.path.join(dirPath, 'cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.vcf')

# Read VCF file
header_vcf = pandas.read_csv(vcfPath, sep = '\t', header = None, index_col = False, dtype = object, nrows = 3473)
vcf = pandas.read_csv(vcfPath, sep = '\t', index_col = False, dtype = object, skiprows = 3474)
vcf = vcf.rename(columns = {'#CHROM': 'CHROM'})
print('Lenght of original vcf')
print(len(vcf))


# Obtain annotation from INFO column
infoCol = vcf['INFO']

## Annotation is in each row after CSQ=
annotation = infoCol.str.extract('CSQ=(.*)')
annotation.columns = ['annotationCol']


# Split annotation by ',' and separate in different rows
annotationSplitComma = annotation['annotationCol'].str.split(',').apply(pandas.Series, 1).stack()
annotationSplitComma.index = annotationSplitComma.index.droplevel(-1)
annotationSplitComma.name = 'annotationCol'

# Create the annotated DF
## Drop not desireed columns
annotatedDF = vcf.drop(['QUAL', 'FILTER', 'INFO', 'FORMAT'], axis = 1)
# Joint the original DF with the DF with the annotation not split by commas
annotatedDFsplitComma = annotatedDF.join(annotationSplitComma)

annotatedDFsplitComma = annotatedDFsplitComma.reset_index() # Reset index to avoid errors
annotatedDFsplitComma = annotatedDFsplitComma.drop(['index'], axis = 1)
print('len annotatedDFsplitComma completed')
print(len(annotatedDFsplitComma))

# Split annotation column by '|'
annotationColSplit = pandas.DataFrame(annotatedDFsplitComma['annotationCol'])
## Count the times that '|' appears, which is the delimeter to split each annotation
counter = Counter(annotationColSplit.iloc[0,0])
counter = counter['|']
print('counter |')
print(counter)
## Split each annotation in columns
annotationColSplit = annotationColSplit['annotationCol'].str.split('|', counter, expand = True)


## Add colnames
### Obtain colnames from header, the names are between Format: and "> in the row that starts with ##INFO=<ID=CSQ
headerAnnotation = header_vcf[header_vcf[0].str.contains('##INFO=<ID=CSQ')]
headerAnnotation = headerAnnotation[0].str.extract('Format: (.*)">')
headerAnnotationSplit = headerAnnotation[0].str.split('|', counter, expand = True)
headerAnnotationSplit = headerAnnotationSplit.values.tolist()[0]

annotationColSplit.columns = headerAnnotationSplit

# Save those to order the DF
colsOrder = annotatedDF.columns.values.tolist()[0:5]
colsOrder.extend(annotationColSplit.columns.values.tolist())
colsOrder.extend(annotatedDF.columns.values.tolist()[5:])
print('Columns order to save')
print(colsOrder)

# Join with rest DF
annotatedDF = pandas.concat([annotatedDFsplitComma, annotationColSplit], axis = 1)

## Reorder columns
annotatedDF = annotatedDF[colsOrder]
print('len final df')
print(len(annotatedDF))

# Write files
nameOutAnnotation = os.path.splitext(vcfPath)[0]
nameOutAnnotation = '.'.join([nameOutAnnotation, 'withAmpersanInCSVS', 'tsv'])
print(nameOutAnnotation)
annotatedDF.to_csv(nameOutAnnotation, header = True, index = None, sep = '\t')

# In CSVS columns, change cells with '&' to NaN
annotatedDF_noAmp = annotatedDF

## Select those columns that contains CSVS in the colname
colsCSVS = [c for c in annotatedDF_noAmp.columns.values.tolist() if 'CSVS' in c]
print(colsCSVS)

for c in colsCSVS:
    annotatedDF_noAmp.loc[annotatedDF_noAmp[c].str.contains('&', na = False), c] = numpy.nan


# Write files without &
nameOutAnnotation_noAmp = os.path.splitext(vcfPath)[0]
nameOutAnnotation_noAmp = '.'.join([nameOutAnnotation_noAmp, 'tsv'])
print(nameOutAnnotation_noAmp)
annotatedDF_noAmp.to_csv(nameOutAnnotation_noAmp, header = True, index = None, sep = '\t')
