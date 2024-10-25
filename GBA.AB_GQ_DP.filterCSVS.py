import pandas
import os
import math
import numpy
from scipy import stats
from pathlib import Path

pandas.set_option('display.max_rows', 200, 'display.max_columns', 1000) # change it to see more or less rows and/or columns

# Define function to perform the GBA
def GBA_function(dirPath, inputName, outPath, inputAnalysesPath):
    # Set paths to save each file
    outPath_ast = os.path.join(dirPath, 'withAsterisc')
    outPath_complete = os.path.join(outPath, 'completeGBA')
    outPath_forDB = os.path.join(outPath, 'for_each_DB')
    outPath_forDB_variants = os.path.join(outPath, 'for_each_DB', 'variants')
    outPath_enriched = os.path.join(outPath, 'at_least_one_DB_enriched')
    outPath_enriched_variants = os.path.join(outPath, 'at_least_one_DB_enriched', 'variants')
    outPath_novel = os.path.join(outPath, 'novelSNV')
    Path(outPath_ast).mkdir(parents = True, exist_ok = True)
    Path(outPath_complete).mkdir(parents = True, exist_ok = True)
    Path(outPath_forDB).mkdir(parents = True, exist_ok = True)
    Path(outPath_forDB_variants).mkdir(parents = True, exist_ok = True)
    Path(outPath_enriched).mkdir(parents = True, exist_ok = True)
    Path(outPath_enriched_variants).mkdir(parents = True, exist_ok = True)
    Path(outPath_novel).mkdir(parents = True, exist_ok = True)


    # Read input file
    inputTSV = pandas.read_csv(os.path.join(dirPath, inputName), sep = '\t', index_col = False, dtype = object)
    print('Number variants in input: ' + str(len(inputTSV)))

    # Filter variants if ALT == *
    inputTSV_ast = inputTSV[inputTSV['ALT'] == '*']
    inputTSV = inputTSV[inputTSV['ALT'] != '*']

    # Write files
    outName = os.path.join(outPath_ast, os.path.splitext(inputName)[0])
    outName_ast = '.'.join([outName, 'withAsterisc', 'tsv'])
    inputTSV_ast.to_csv(outName_ast, header = True, index = None, sep = '\t', float_format = '%f')

    # Transform some columns to float
    colsToFloat = ['CADD_PHRED', 'CADD_RAW', 'gnomADg_AF_nfe', 'gnomADg_AF', 'CSVS_AF', 'AF']
    for col in colsToFloat:
        inputTSV[col] = inputTSV[col].astype(float)

    dtypes = dict(inputTSV.dtypes) # Check the types of each column


    # Filter variants without symbol
    inputTSV_SYMBOL = inputTSV[inputTSV['SYMBOL'].notnull()]
    print('Variants with gene name: ' + str(len(inputTSV_SYMBOL)))


    # Prepare data to filter by AF
    ## SNVs without data for gnomADg and CSVS should have it to perform the GBA with values for these variants
    ## Variant without data will be: AC = 0, AF = 0, AN = mean(AN in gene)
    ### Create table with meanAN for each genes (all SNVs, before filters)
    #### gnomAD
    cols_meanAN_gnomAD = ['SYMBOL', 'gnomADg_AN_nfe', 'gnomADg_AN']
    meanAN_DF_gnomAD = inputTSV[cols_meanAN_gnomAD] # DF with less columns
    meanAN_DF_gnomAD = meanAN_DF_gnomAD[meanAN_DF_gnomAD['gnomADg_AN_nfe'].notnull()] # Variants with value
    print('Variants with value in gnomAD: ' + str(len(meanAN_DF_gnomAD)))
    colsToNumeric_gnomAD = ['gnomADg_AN_nfe', 'gnomADg_AN'] # Transform columns to numeric
    for col in colsToNumeric_gnomAD:
        meanAN_DF_gnomAD[col] = meanAN_DF_gnomAD[col].astype('int32')
    meanAN_DF_gnomAD = meanAN_DF_gnomAD.groupby(by = ['SYMBOL']).mean() # Calculate the mean for each gene
    meanAN_DF_gnomAD = meanAN_DF_gnomAD.round(0).astype(int) # Number without decimals
    meanAN_DF_gnomAD = meanAN_DF_gnomAD.reset_index() # Reset index to avoid errors
    print('Genes with value in gnomAD: ' + str(len(meanAN_DF_gnomAD)))

    #### CSVS
    cols_meanAN_CSVS = ['SYMBOL', 'CSVS_AN']
    meanAN_DF_CSVS = inputTSV[cols_meanAN_CSVS]
    meanAN_DF_CSVS = meanAN_DF_CSVS[meanAN_DF_CSVS['CSVS_AN'].notnull()]
    print('Variants with value in CSVS: ' + str(len(meanAN_DF_CSVS)))
    colsToNumeric_CSVS = ['CSVS_AN']
    for col in colsToNumeric_CSVS:
        meanAN_DF_CSVS[col] = meanAN_DF_CSVS[col].astype('int32')
    meanAN_DF_CSVS = meanAN_DF_CSVS.groupby(by = ['SYMBOL']).mean()
    meanAN_DF_CSVS = meanAN_DF_CSVS.round(0).astype(int)
    meanAN_DF_CSVS = meanAN_DF_CSVS.reset_index()
    print('Genes with value in CSVS: ' + str(len(meanAN_DF_CSVS)))

    ### Merge gnomAD and CSVS
    meanAN_DF = pandas.merge(meanAN_DF_gnomAD, meanAN_DF_CSVS, how = 'left', left_on = 'SYMBOL', right_on = 'SYMBOL')

    #### Genes without value in both databases
    symbol_diff = list(set(inputTSV.SYMBOL.unique().tolist()).difference(set(meanAN_DF.SYMBOL.unique().tolist())))
    symbol_diffDF = pandas.DataFrame(columns = meanAN_DF.columns.values.tolist(), index = range(len(symbol_diff)))
    symbol_diffDF['SYMBOL'] = symbol_diff
    meanAN_DF = pandas.concat([meanAN_DF, symbol_diffDF], ignore_index=True)

    ### If there is a gene without any value use as AN the mean of all the genes for the same database
    for c in meanAN_DF.columns[1:]:
        for r in meanAN_DF.index:
            if math.isnan(meanAN_DF.loc[r,c]):
                # print('There is a nan in the column: ' + str(c) + ' and in the row: ' + str(r))
                meanAN_col = numpy.mean(meanAN_DF[c])
                meanAN_DF.loc[r,c] = meanAN_col.round(0).astype(int)
    colsToNumeric_meanAN_DF = ['gnomADg_AN_nfe', 'gnomADg_AN', 'CSVS_AN'] # Transform columns to numeric
    for col in colsToNumeric_meanAN_DF:
        meanAN_DF[col] = meanAN_DF[col].astype('int32')

    ## Add values to DF
    inputTSV_AF = inputTSV_SYMBOL
    for r in inputTSV_AF.index:
        if pandas.isnull(inputTSV_AF['gnomADg'][r]) == True:
            # AF & AC
            inputTSV_AF.at[r, 'gnomADg_AF_nfe'] = 0
            inputTSV_AF.at[r, 'gnomADg_AC_nfe'] = 0
            inputTSV_AF.at[r, 'gnomADg_AF'] = 0
            inputTSV_AF.at[r, 'gnomADg_AC'] = 0
            # AN
            rGene = inputTSV_AF['SYMBOL'][r]
            ANGene_nfe = meanAN_DF.loc[meanAN_DF['SYMBOL'] == rGene, 'gnomADg_AN_nfe'].iloc[0]
            ANGene = meanAN_DF.loc[meanAN_DF['SYMBOL'] == rGene, 'gnomADg_AN'].iloc[0]
            inputTSV_AF.at[r, 'gnomADg_AN_nfe'] = ANGene_nfe
            inputTSV_AF.at[r, 'gnomADg_AN'] = ANGene
        if pandas.isnull(inputTSV_AF['CSVS'][r]) == True:
            # AF & AC
            inputTSV_AF.at[r, 'CSVS_AF'] = 0
            inputTSV_AF.at[r, 'CSVS_AC'] = 0
            # AN
            rGene = inputTSV_AF['SYMBOL'][r]
            ANGene_CSVS = meanAN_DF.loc[meanAN_DF['SYMBOL'] == rGene, 'CSVS_AN'].iloc[0]
            inputTSV_AF.at[r, 'CSVS_AN'] = ANGene_CSVS


    # Show to the user the options to filter by consequence and impact
    print('Consequences in input: ' + ', '.join(inputTSV_AF.Consequence.unique()))
    print('Impacts in input: ' + ', '.join(inputTSV_AF.IMPACT.unique()))


    # Set the filters for each analysis
    inputAnalyses = pandas.read_csv(inputAnalysesPath, sep = '\t', index_col = False, dtype = object)
    print('Number of analyses: ' + str(len(inputAnalyses)))

    ## Transform some columns to float
    colsToFloat = ['CADD', 'AF_nfe', 'AF_global', 'AF_CSVS']
    for col in colsToFloat:
        inputAnalyses[col] = inputAnalyses[col].astype(float)

    dtypes = dict(inputAnalyses.dtypes) # Check the types of each column

    ## Create lists from each column
    snp_indel = inputAnalyses['snp_indel'].values.tolist()
    consequence = inputAnalyses['consequence'].values.tolist()
    impact = inputAnalyses['impact'].values.tolist()
    LOFTEE = inputAnalyses['LOFTEE'].values.tolist()
    CADD = inputAnalyses['CADD'].values.tolist()
    AF_nfe = inputAnalyses['AF_nfe'].values.tolist()
    AF_global = inputAnalyses['AF_global'].values.tolist()
    AF_CSVS = inputAnalyses['AF_CSVS'].values.tolist()


    # Create dataframe to save summary data
    columnsSummary = inputAnalyses.columns.values.tolist()
    columnsSummary = columnsSummary + ['n_genes_consequence','n_variants_consequence','n_genes_impact','n_variants_impact','n_genes_LOFTEE','n_variants_LOFTEE','n_genes_CADD','n_variants_CADD','n_genes_AF','n_variants_AF','n_genes_nfe','n_genes_global','n_genes_nfe_global','n_genes_CSVS','n_genes_nfe_global_CSVS','n_variants_nfe','n_variants_global','n_variants_nfe_global','n_variants_CSVS','n_variants_nfe_global_CSVS','n_genes_oneDB','n_variants_oneDB']
    summaryData = inputAnalyses.reindex(columns=columnsSummary)


    # Perform GBA for each analysis
    for i in range(len(consequence)):
        ## Ask inputs to filter
        snp_indel_filter = snp_indel[i]
        consequence_filter = consequence[i]
        impact_filter = impact[i]
        LOFTEE_filter = LOFTEE[i]
        CADD_filter = CADD[i]
        gnomADg_AF_nfe_filter = AF_nfe[i]
        gnomADg_AF_filter = AF_global[i]
        CSVS_AF_filter = AF_CSVS[i]

        print('Analysis:')
        print(inputAnalyses.loc[i])
    
        # Filter by SNPs/indels
        nt = ['A', 'T', 'C', 'G']
        if snp_indel_filter == 'snp':
            inputTSV_nt_AF = inputTSV_AF[inputTSV_AF['REF'].isin(nt) & inputTSV_AF['ALT'].isin(nt)]
            inputTSV_nt_AF.REF.unique()
            inputTSV_nt_AF.ALT.unique()
        elif snp_indel_filter == 'indel':
            inputTSV_nt_AF = inputTSV_AF[~inputTSV_AF['REF'].isin(nt) | ~inputTSV_AF['ALT'].isin(nt)]
        elif snp_indel_filter == 'snp|indel':
            inputTSV_nt_AF = inputTSV_AF
        else:
            print('Bad snp/indel filter introduced, the options are: "snp", "indel" or "snp|indel"')

        print('Variants after filtering by ' + snp_indel_filter + ': ' + str(len(inputTSV_nt_AF)))


        # Filter by variant type (Consequence)
        if consequence_filter == 'all':
            consequenceDF = inputTSV_nt_AF
        else:
            consequenceDF = inputTSV_nt_AF[inputTSV_nt_AF['Consequence'].str.contains(consequence_filter)]

        print('Genes with variants after filtering by the consequence(s) ' + consequence_filter + ' : ' + str(len(consequenceDF['SYMBOL'].unique())))
        print('Variants after filtering by the consequence(s) ' + consequence_filter + ' : ' + str(len(consequenceDF)))
        summaryData.at[i,'n_genes_consequence'] = len(consequenceDF['SYMBOL'].unique())
        summaryData.at[i,'n_variants_consequence'] = len(consequenceDF)


        # Filter by impact
        if impact_filter == 'all':
            consequence_impactDF = consequenceDF
        else:
            consequence_impactDF = consequenceDF[consequenceDF['IMPACT'].str.contains(impact_filter)]

        print('Genes with variants after filtering by the impact(s) ' + impact_filter + ' : ' + str(len(consequence_impactDF['SYMBOL'].unique())))
        print('Variants after filtering by the impact(s) ' + impact_filter + ' : ' + str(len(consequence_impactDF)))
        summaryData.at[i,'n_genes_impact'] = len(consequence_impactDF['SYMBOL'].unique())
        summaryData.at[i,'n_variants_impact'] = len(consequence_impactDF)


        # Filter by LOFTEE
        if pandas.isna(LOFTEE_filter):
            consequence_impact_LOFTEEDF = consequence_impactDF
            LOFTEE_filter_name = ''
            print('LOFTEE filter is NaN')
        else:
            consequence_impact_LOFTEEDF = consequence_impactDF[((consequence_impactDF['LoF'] == LOFTEE_filter) & (consequence_impactDF['IMPACT'] == 'HIGH')) | (consequence_impactDF['IMPACT'] != 'HIGH')]
            print('Genes with variants after filtering by the LOFTEE ' + LOFTEE_filter + ' : ' + str(len(consequence_impact_LOFTEEDF['SYMBOL'].unique())))
            print('Variants after filtering by the LOFTEE ' + LOFTEE_filter + ' : ' + str(len(consequence_impact_LOFTEEDF)))
            summaryData.at[i,'n_genes_LOFTEE'] = len(consequence_impact_LOFTEEDF['SYMBOL'].unique())
            summaryData.at[i,'n_variants_LOFTEE'] = len(consequence_impact_LOFTEEDF)
            LOFTEE_filter_name = LOFTEE_filter


        # Filter by CADD
        if pandas.isna(CADD_filter):
            consequence_impact_LOFTEE_CADDDF = consequence_impact_LOFTEEDF
            CADD_filter_name = ''
            print('CADD filter is NaN')
        else:
            consequence_impact_LOFTEE_CADDDF = consequence_impact_LOFTEEDF[((consequence_impact_LOFTEEDF['CADD_PHRED'] >= CADD_filter) & (consequence_impact_LOFTEEDF['IMPACT'] == 'MODERATE')) | (consequence_impact_LOFTEEDF['IMPACT'] != 'MODERATE')]
            print('Genes with variants after filtering by the CADD ' + str(CADD_filter) + ' : ' + str(len(consequence_impact_LOFTEE_CADDDF['SYMBOL'].unique())))
            print('Variants after filtering by the CADD ' + str(CADD_filter) + ' : ' + str(len(consequence_impact_LOFTEE_CADDDF)))
            summaryData.at[i,'n_genes_CADD'] = len(consequence_impact_LOFTEE_CADDDF['SYMBOL'].unique())
            summaryData.at[i,'n_variants_CADD'] = len(consequence_impact_LOFTEE_CADDDF)
            CADD_filter_name = CADD_filter

    
        # Filter by AF of gnomAD nfe, gnomAD and CSVS
        consequence_impact_LOFTEE_CADD_AFDF = consequence_impact_LOFTEE_CADDDF[(consequence_impact_LOFTEE_CADDDF['gnomADg_AF_nfe'] < gnomADg_AF_nfe_filter) & (consequence_impact_LOFTEE_CADDDF['gnomADg_AF'] < gnomADg_AF_filter) & (consequence_impact_LOFTEE_CADDDF['CSVS_AF'] < CSVS_AF_filter)]
        print('Genes with variants after filtering by the AF ' + str(gnomADg_AF_nfe_filter) + ', ' + str(gnomADg_AF_filter) + ', ' +  str(CSVS_AF_filter) + ' : ' + str(len(consequence_impact_LOFTEE_CADD_AFDF['SYMBOL'].unique())))
        print('Variants after filtering by the AF ' + str(gnomADg_AF_nfe_filter) + ', ' + str(gnomADg_AF_filter) + ', ' +  str(CSVS_AF_filter) + ' : ' + str(len(consequence_impact_LOFTEE_CADD_AFDF)))
        summaryData.at[i,'n_genes_AF'] = len(consequence_impact_LOFTEE_CADD_AFDF['SYMBOL'].unique())
        summaryData.at[i,'n_variants_AF'] = len(consequence_impact_LOFTEE_CADD_AFDF)


        # Check if exist variants
        if len(consequence_impact_LOFTEE_CADD_AFDF) == 0:
            continue
        else:
            print('After these filters there are ' + str(len(consequence_impact_LOFTEE_CADD_AFDF)) + ' variants to continue with the GBA\n')



        # Create DF as GBA input
        ## The necessary columns are the gene and the AC and AN for the cases and each database
        colsGBA = ['SYMBOL', 'AC', 'AN', 'gnomADg_AC_nfe', 'gnomADg_AN_nfe', 'gnomADg_AC', 'gnomADg_AN', 'CSVS_AC', 'CSVS_AN'] # AN are the total number of alleles
        inputGBA_SNV = consequence_impact_LOFTEE_CADD_AFDF[colsGBA]

        ## Rename cases colnames
        colsRenameIndex = [1,2]
        namesColsNew = ['AC_cases', 'AN_cases']
        namesColsOld = inputGBA_SNV.columns[colsRenameIndex]
        inputGBA_SNV.rename(columns = dict(zip(namesColsOld, namesColsNew)), inplace = True)

        ## Change to integer, except SYMBOL column
        for col in inputGBA_SNV.columns[1:]:
            inputGBA_SNV[col] = inputGBA_SNV[col].astype(int)

        dtypes = dict(inputGBA_SNV.dtypes) # Check the types of each column


        ## Calculate WT (wild type): WT = total alleles - allele count (variant)
        inputGBA_SNV['WT_cases'] = inputGBA_SNV['AN_cases'] - inputGBA_SNV['AC_cases']
        inputGBA_SNV['WT_gnomADg_nfe'] = inputGBA_SNV['gnomADg_AN_nfe'] - inputGBA_SNV['gnomADg_AC_nfe']
        inputGBA_SNV['WT_gnomADg'] = inputGBA_SNV['gnomADg_AN'] - inputGBA_SNV['gnomADg_AC']
        inputGBA_SNV['WT_CSVS'] = inputGBA_SNV['CSVS_AN'] - inputGBA_SNV['CSVS_AC']

        ## Remove columns with AN
        inputGBA_SNV = inputGBA_SNV[inputGBA_SNV.columns.drop(list(inputGBA_SNV.filter(regex = 'AN')))]

        ## Calculate the sum for each column grouping by gene name
        inputGBA = inputGBA_SNV.groupby(by = ['SYMBOL']).sum()
        inputGBA = inputGBA.reset_index()
        print('The number of genes for the GBA is: ' + str(len(inputGBA)))


        # Extract genes with all SNV novel
        ## Is not possible to perform a gene burden with values = 0
        ## Study separatly
        novelSNV_genes = inputGBA[(inputGBA['gnomADg_AC_nfe'] == 0) | (inputGBA['gnomADg_AC'] == 0) | (inputGBA['CSVS_AC'] == 0)]


        # Odd Ratio
        inputGBA_OR = inputGBA
        inputGBA_OR['OR_gnomADg_nfe'] = (inputGBA_OR['AC_cases']*inputGBA_OR['WT_gnomADg_nfe'])/(inputGBA_OR['WT_cases']*inputGBA_OR['gnomADg_AC_nfe'])
        inputGBA_OR['OR_gnomADg'] = (inputGBA_OR['AC_cases']*inputGBA_OR['WT_gnomADg'])/(inputGBA_OR['WT_cases']*inputGBA_OR['gnomADg_AC'])
        inputGBA_OR['OR_CSVS'] = (inputGBA_OR['AC_cases']*inputGBA_OR['WT_CSVS'])/(inputGBA_OR['WT_cases']*inputGBA_OR['CSVS_AC'])

        # Standard error
        ## Calculate the summary
        inputGBA_OR_SE = inputGBA_OR
        sumList_nfe = ((1/inputGBA_OR_SE['AC_cases']) + (1/inputGBA_OR_SE['gnomADg_AC_nfe']) + (1/inputGBA_OR_SE['WT_cases']) + (1/inputGBA_OR_SE['WT_gnomADg_nfe'])).tolist()
        sumList_gnomAD = ((1/inputGBA_OR_SE['AC_cases']) + (1/inputGBA_OR_SE['gnomADg_AC']) + (1/inputGBA_OR_SE['WT_cases']) + (1/inputGBA_OR_SE['WT_gnomADg'])).tolist()
        sumList_CSVS = ((1/inputGBA_OR_SE['AC_cases']) + (1/inputGBA_OR_SE['CSVS_AC']) + (1/inputGBA_OR_SE['WT_cases']) + (1/inputGBA_OR_SE['WT_CSVS'])).tolist()
        ## Perform the sqrt
        SElist_nfe = []
        for sum in sumList_nfe:
            SE = math.sqrt(sum)
            SElist_nfe.append(SE)

        SElist_gnomAD = []
        for sum in sumList_gnomAD:
            SE = math.sqrt(sum)
            SElist_gnomAD.append(SE)

        SElist_CSVS = []
        for sum in sumList_CSVS:
            SE = math.sqrt(sum)
            SElist_CSVS.append(SE)


        inputGBA_OR_SE['SE_gnomADg_nfe'] = SElist_nfe
        inputGBA_OR_SE['SE_gnomADg'] = SElist_gnomAD
        inputGBA_OR_SE['SE_CSVS'] = SElist_CSVS
        # inputGBA['SE_gnomADg_nfe'] = math.sqrt((1/inputGBA['AC_cases']) + (1/inputGBA['gnomADg_AC_nfe']) + (1/inputGBA['WT_cases']) + (1/inputGBA['WT_gnomADg_nfe'])) --> doesn't work


        # Z-Score
        inputGBA_OR_SE_Z = inputGBA_OR_SE

        inputGBA_OR_SE_Z['ZScore_gnomADg_nfe'] = numpy.log(inputGBA_OR_SE_Z['OR_gnomADg_nfe'])/inputGBA_OR_SE_Z['SE_gnomADg_nfe']
        inputGBA_OR_SE_Z['ZScore_gnomADg'] = numpy.log(inputGBA_OR_SE_Z['OR_gnomADg'])/inputGBA_OR_SE_Z['SE_gnomADg']
        inputGBA_OR_SE_Z['ZScore_CSVS'] = numpy.log(inputGBA_OR_SE_Z['OR_CSVS'])/inputGBA_OR_SE_Z['SE_CSVS']


        # p-value
        inputGBA_OR_SE_Z_pv = inputGBA_OR_SE_Z

        inputGBA_OR_SE_Z_pv['pvalue_gnomADg_nfe'] = stats.norm.sf(abs(inputGBA_OR_SE_Z['ZScore_gnomADg_nfe']))*2 # using CCDF
        inputGBA_OR_SE_Z_pv['pvalue_gnomADg'] = stats.norm.sf(abs(inputGBA_OR_SE_Z['ZScore_gnomADg']))*2
        inputGBA_OR_SE_Z_pv['pvalue_CSVS'] = stats.norm.sf(abs(inputGBA_OR_SE_Z['ZScore_CSVS']))*2
        ### (1 - stats.norm.cdf(abs(inputGBA_OR_SE_Z['ZScore_gnomADg_nfe'])))*2 # using CDF --> same: 1 - CDF = CCDF

        # FDR - number of genes
        inputGBA_OR_SE_Z_pv['FDR_gnomADg_nfe'] = inputGBA_OR_SE_Z_pv['pvalue_gnomADg_nfe'] * len(inputGBA_OR_SE_Z_pv[inputGBA_OR_SE_Z_pv['gnomADg_AC_nfe'] != 0]) # number of genes in the analysis
        inputGBA_OR_SE_Z_pv['FDR_gnomADg'] = inputGBA_OR_SE_Z_pv['pvalue_gnomADg'] * len(inputGBA_OR_SE_Z_pv[inputGBA_OR_SE_Z_pv['gnomADg_AC'] != 0])
        inputGBA_OR_SE_Z_pv['FDR_CSVS'] = inputGBA_OR_SE_Z_pv['pvalue_CSVS'] * len(inputGBA_OR_SE_Z_pv[inputGBA_OR_SE_Z_pv['CSVS_AC'] != 0])

        # Confidence interval
        inputGBA_OR_SE_Z_pv_CI = inputGBA_OR_SE_Z_pv
        inputGBA_OR_SE_Z_pv_CI['lowCI_gnomADg_nfe'] = numpy.exp(numpy.log(inputGBA_OR_SE_Z_pv['OR_gnomADg_nfe']) - 1.96*inputGBA_OR_SE_Z_pv['SE_gnomADg_nfe'])
        inputGBA_OR_SE_Z_pv_CI['highCI_gnomADg_nfe'] = numpy.exp(numpy.log(inputGBA_OR_SE_Z_pv['OR_gnomADg_nfe']) + 1.96*inputGBA_OR_SE_Z_pv['SE_gnomADg_nfe'])
        inputGBA_OR_SE_Z_pv_CI['lowCI_gnomADg'] = numpy.exp(numpy.log(inputGBA_OR_SE_Z_pv['OR_gnomADg']) - 1.96*inputGBA_OR_SE_Z_pv['SE_gnomADg'])
        inputGBA_OR_SE_Z_pv_CI['highCI_gnomADg'] = numpy.exp(numpy.log(inputGBA_OR_SE_Z_pv['OR_gnomADg']) + 1.96*inputGBA_OR_SE_Z_pv['SE_gnomADg'])
        inputGBA_OR_SE_Z_pv_CI['lowCI_CSVS'] = numpy.exp(numpy.log(inputGBA_OR_SE_Z_pv['OR_CSVS']) - 1.96*inputGBA_OR_SE_Z_pv['SE_CSVS'])
        inputGBA_OR_SE_Z_pv_CI['highCI_CSVS'] = numpy.exp(numpy.log(inputGBA_OR_SE_Z_pv['OR_CSVS']) + 1.96*inputGBA_OR_SE_Z_pv['SE_CSVS'])

        # Etiological fraction
        inputGBA_OR_SE_Z_pv_CI_EF = inputGBA_OR_SE_Z_pv_CI
        inputGBA_OR_SE_Z_pv_CI_EF['EF_gnomADg_nfe'] = (inputGBA_OR_SE_Z_pv_CI_EF['OR_gnomADg_nfe'] - 1) / inputGBA_OR_SE_Z_pv_CI_EF['OR_gnomADg_nfe']
        inputGBA_OR_SE_Z_pv_CI_EF['EF_gnomADg'] = (inputGBA_OR_SE_Z_pv_CI_EF['OR_gnomADg'] - 1) / inputGBA_OR_SE_Z_pv_CI_EF['OR_gnomADg']
        inputGBA_OR_SE_Z_pv_CI_EF['EF_CSVS'] = (inputGBA_OR_SE_Z_pv_CI_EF['OR_CSVS'] - 1) / inputGBA_OR_SE_Z_pv_CI_EF['OR_CSVS']


        # Add variants and samples to the result of the GBA
        GBA = inputGBA_OR_SE_Z_pv_CI_EF
        for r in GBA.index:
            # print(r)
            gene = GBA.iloc[r,0]
            variants_gene = consequence_impact_LOFTEE_CADD_AFDF[consequence_impact_LOFTEE_CADD_AFDF['SYMBOL'] == gene]
            ## Change NaN in samples_het and samples_hom columns to avoid errors
            variants_gene['samples_het'] = variants_gene['samples_het'].replace(numpy.NaN, '')
            variants_gene['samples_hom'] = variants_gene['samples_hom'].replace(numpy.NaN, '')
            ## Count number of variants in the gene
            nVars = len(variants_gene)
            ## Number and name of heterocygous samples
            samps_het = '|'.join(variants_gene['samples_het'].values.tolist())
            samps_het = samps_het.split('|')
            samps_het = list(set(samps_het))
            samps_het = list(filter(None, samps_het))
            n_samps_het = len(samps_het)
            samps_het = '|'.join(samps_het)
            ## Number and name of homorocygous samples
            samps_hom = '|'.join(variants_gene['samples_hom'].values.tolist())
            samps_hom = samps_hom.split('|')
            samps_hom = list(set(samps_hom))
            samps_hom = list(filter(None, samps_hom))
            n_samps_hom = len(samps_hom)
            samps_hom = '|'.join(samps_hom)

            # Add values in each column for the gene
            GBA.at[r, 'n_variants'] = nVars
            GBA.at[r, 'n_individuals'] = n_samps_het + n_samps_hom
            GBA.at[r, 'n_het'] = n_samps_het
            GBA.at[r, 'n_hom'] = n_samps_hom
            GBA.at[r, 'samples_het'] = samps_het
            GBA.at[r, 'samples_hom'] = samps_hom

        ## Convert these columns to integer
        GBA[['n_variants','n_individuals','n_het','n_hom']] = GBA[['n_variants','n_individuals','n_het','n_hom']].astype(int)

        # Reorder columns
        colsOrderFinal = ['SYMBOL', 'n_variants','n_individuals','n_het','n_hom', 'AC_cases', 'WT_cases', 'gnomADg_AC_nfe', 'WT_gnomADg_nfe', 'OR_gnomADg_nfe', 'SE_gnomADg_nfe', 'ZScore_gnomADg_nfe', 'pvalue_gnomADg_nfe', 'FDR_gnomADg_nfe', 'lowCI_gnomADg_nfe', 'highCI_gnomADg_nfe', 'EF_gnomADg_nfe', 'gnomADg_AC', 'WT_gnomADg',	'OR_gnomADg', 'SE_gnomADg',	'ZScore_gnomADg', 'pvalue_gnomADg', 'FDR_gnomADg', 'lowCI_gnomADg', 'highCI_gnomADg', 'EF_gnomADg', 'CSVS_AC', 'WT_CSVS',	'OR_CSVS', 'SE_CSVS',	'ZScore_CSVS', 'pvalue_CSVS', 'FDR_CSVS', 'lowCI_CSVS', 'highCI_CSVS', 'EF_CSVS', 'samples_het', 'samples_hom']
        GBA = GBA[colsOrderFinal]


        # Filter by FDR < 0.05 in each database and combining them
        ## In CSVS, NaN values are allowed, because the database is small
        GBA_nfe = GBA[GBA['FDR_gnomADg_nfe'] < 0.05] # Enriched when comparing with NFE
        GBA_global = GBA[GBA['FDR_gnomADg'] < 0.05] # Enriched when comparing with global
        GBA_nfe_global = GBA[(GBA['FDR_gnomADg_nfe'] < 0.05) & (GBA['FDR_gnomADg'] < 0.05)] # Enriched when comparing with NFE and global
        GBA_CSVS = GBA[(GBA['FDR_CSVS'] < 0.05) | (GBA['FDR_CSVS'].isnull())] # Enriched when comparing with CSVS or without data in CSVS
        GBA_all = GBA[(GBA['FDR_gnomADg_nfe'] < 0.05) & (GBA['FDR_gnomADg'] < 0.05) & ((GBA['FDR_CSVS'] < 0.05) | (GBA['FDR_CSVS'].isnull()))] # Enriched when comparing with NFE, global and CSVS or without data in CSVS
        GBA_enriched = GBA[(GBA['FDR_gnomADg_nfe'] < 0.05) | (GBA['FDR_gnomADg'] < 0.05) | ((GBA['FDR_CSVS'] < 0.05) | (GBA['FDR_CSVS'].isnull()))] # Enriched at least when comparing with one of: NFE, global, or CSVS or without data in CSVS


        # Print the results of the GBA
        print('The GBA is done. Filtering by FDR < 0.05 for the following databases, the next number of genes were enriched:')
        print('gnomAD NFE: ' + str(len(GBA_nfe)) + ' genes')
        print('gnomAD global population: ' + str(len(GBA_global)) + ' genes')
        print('gnomAD NFE + gnomAD global population: ' + str(len(GBA_nfe_global)) + ' genes')
        print('CSVS: ' + str(len(GBA_CSVS)) + ' genes')
        print('gnomAD NFE + gnomAD global population + CSVS: ' + str(len(GBA_all)) + ' genes')
        print('at least one of them: ' + str(len(GBA_enriched)) + ' genes')

        summaryData.at[i,'n_genes_nfe'] = len(GBA_nfe)
        summaryData.at[i,'n_genes_global'] = len(GBA_global)
        summaryData.at[i,'n_genes_nfe_global'] = len(GBA_nfe_global)
        summaryData.at[i,'n_genes_CSVS'] = len(GBA_CSVS)
        summaryData.at[i,'n_genes_nfe_global_CSVS'] = len(GBA_all)
        summaryData.at[i,'n_genes_oneDB'] = len(GBA_enriched)


        # Variants in enriched genes
        GBA_nfe_variants = consequence_impact_LOFTEE_CADD_AFDF[consequence_impact_LOFTEE_CADD_AFDF['SYMBOL'].isin(GBA_nfe['SYMBOL'].values.tolist())]
        GBA_global_variants = consequence_impact_LOFTEE_CADD_AFDF[consequence_impact_LOFTEE_CADD_AFDF['SYMBOL'].isin(GBA_global['SYMBOL'].values.tolist())]
        GBA_nfe_global_variants = consequence_impact_LOFTEE_CADD_AFDF[consequence_impact_LOFTEE_CADD_AFDF['SYMBOL'].isin(GBA_nfe_global['SYMBOL'].values.tolist())]
        GBA_CSVS_variants = consequence_impact_LOFTEE_CADD_AFDF[consequence_impact_LOFTEE_CADD_AFDF['SYMBOL'].isin(GBA_CSVS['SYMBOL'].values.tolist())]
        GBA_all_variants = consequence_impact_LOFTEE_CADD_AFDF[consequence_impact_LOFTEE_CADD_AFDF['SYMBOL'].isin(GBA_all['SYMBOL'].values.tolist())]
        GBA_enriched_variants = consequence_impact_LOFTEE_CADD_AFDF[consequence_impact_LOFTEE_CADD_AFDF['SYMBOL'].isin(GBA_enriched['SYMBOL'].values.tolist())]

        print('Number of variants in enriched genes in gnomAD NFE: ', str(len(GBA_nfe_variants)))
        print('Number of variants in enriched genes in gnomAD all population: ', str(len(GBA_global_variants)))
        print('Number of variants in enriched genes in gnomAD NFE + gnomAD all population: ', str(len(GBA_nfe_global_variants)))
        print('Number of variants in enriched genes in CSVS: ', str(len(GBA_CSVS_variants)))
        print('Number of variants in enriched genes in gnomAD NFE + gnomAD all population + CSVS: ', str(len(GBA_all_variants)))
        print('Number of variants in enriched genes in at least one DB: ', str(len(GBA_enriched_variants)))

        summaryData.at[i,'n_variants_nfe'] = len(GBA_nfe_variants)
        summaryData.at[i,'n_variants_global'] = len(GBA_global_variants)
        summaryData.at[i,'n_variants_nfe_global'] = len(GBA_nfe_global_variants)
        summaryData.at[i,'n_variants_CSVS'] = len(GBA_CSVS_variants)
        summaryData.at[i,'n_variants_nfe_global_CSVS'] = len(GBA_all_variants)
        summaryData.at[i,'n_variants_oneDB'] = len(GBA_enriched_variants)


        # Write files, using input filters in the file name
        print('Files will be saved in:\n' + outPath_complete + '\n' + outPath_forDB + '\n' + outPath_forDB_variants + '\n' + outPath_enriched + '\n' + outPath_enriched_variants + '\n' + outPath_novel)

        ## Set base name to save
        baseNameSave = '.'.join([os.path.splitext(inputName)[0],'GBA', snp_indel_filter.replace('|', '_'), consequence_filter.replace('|', '_'), impact_filter.replace('|', '_'), LOFTEE_filter_name.replace('|', '_'), str(CADD_filter_name).replace('.', ''), str(gnomADg_AF_nfe_filter).replace('.', ''), str(gnomADg_AF_filter).replace('.', ''), str(CSVS_AF_filter).replace('.', '')])
        print('Base name of files: ' + baseNameSave)

        ## Save files
        outName_GBA = os.path.join(outPath_complete, '.'.join([baseNameSave, 'tsv']))
        GBA.to_csv(outName_GBA, header = True, index = None, sep = '\t', float_format = '%.16f')

        if len(GBA_nfe) != 0:
            outName_GBA_nfe = os.path.join(outPath_forDB, '.'.join([baseNameSave, 'nfe', 'tsv']))
            GBA_nfe.to_csv(outName_GBA_nfe, header = True, index = None, sep = '\t', float_format = '%.16f')
            outName_GBA_nfe_variants = os.path.join(outPath_forDB_variants, '.'.join([baseNameSave, 'nfe', 'variants', 'tsv']))
            GBA_nfe_variants.to_csv(outName_GBA_nfe_variants, header = True, index = None, sep = '\t', float_format = '%.16f')

        if len(GBA_global) != 0:
            outName_GBA_global = os.path.join(outPath_forDB, '.'.join([baseNameSave, 'global', 'tsv']))
            GBA_global.to_csv(outName_GBA_global, header = True, index = None, sep = '\t', float_format = '%.16f')
            outName_GBA_global_variants = os.path.join(outPath_forDB_variants, '.'.join([baseNameSave, 'global', 'variants', 'tsv']))
            GBA_global_variants.to_csv(outName_GBA_global_variants, header = True, index = None, sep = '\t', float_format = '%.16f')

        if len(GBA_nfe_global) != 0:
            outName_GBA_nfe_global = os.path.join(outPath_forDB, '.'.join([baseNameSave, 'nfe', 'global', 'tsv']))
            GBA_nfe_global.to_csv(outName_GBA_nfe_global, header = True, index = None, sep = '\t', float_format = '%.16f')
            outName_GBA_nfe_global_variants = os.path.join(outPath_forDB_variants, '.'.join([baseNameSave, 'nfe', 'global', 'variants', 'tsv']))
            GBA_nfe_global_variants.to_csv(outName_GBA_nfe_global_variants, header = True, index = None, sep = '\t', float_format = '%.16f')

        if len(GBA_CSVS) != 0:
            outName_GBA_CSVS = os.path.join(outPath_forDB, '.'.join([baseNameSave, 'CSVS', 'tsv']))
            GBA_CSVS.to_csv(outName_GBA_CSVS, header = True, index = None, sep = '\t', float_format = '%.16f')
            outName_GBA_CSVS_variants = os.path.join(outPath_forDB_variants, '.'.join([baseNameSave, 'CSVS', 'variants', 'tsv']))
            GBA_CSVS_variants.to_csv(outName_GBA_CSVS_variants, header = True, index = None, sep = '\t', float_format = '%.16f')

        if len(GBA_all) != 0:
            outName_GBA_all = os.path.join(outPath_forDB, '.'.join([baseNameSave, 'nfe', 'global', 'CSVS', 'tsv']))
            GBA_all.to_csv(outName_GBA_all, header = True, index = None, sep = '\t', float_format = '%.16f')
            outName_GBA_all_variants = os.path.join(outPath_forDB_variants, '.'.join([baseNameSave, 'nfe', 'global', 'CSVS', 'variants', 'tsv']))
            GBA_all_variants.to_csv(outName_GBA_all_variants, header = True, index = None, sep = '\t', float_format = '%.16f')

        if len(GBA_enriched) != 0:
            outName_GBA_enriched = os.path.join(outPath_enriched, '.'.join([baseNameSave, 'oneDB', 'tsv']))
            GBA_enriched.to_csv(outName_GBA_enriched, header = True, index = None, sep = '\t', float_format = '%.16f')
            outName_GBA_enriched_variants = os.path.join(outPath_enriched_variants, '.'.join([baseNameSave, 'oneDB', 'variants', 'tsv']))
            GBA_enriched_variants.to_csv(outName_GBA_enriched_variants, header = True, index = None, sep = '\t', float_format = '%.16f')

        if len(novelSNV_genes) != 0:
            outName_GBA_enriched = os.path.join(outPath_novel, '.'.join([baseNameSave, 'novelSNV', 'tsv']))
            novelSNV_genes.to_csv(outName_GBA_enriched, header = True, index = None, sep = '\t', float_format = '%.16f')

        print('All files saved')

    baseNameSaveSummary = '.'.join([os.path.splitext(inputName)[0],'GBA'])
    outName_summary = os.path.join(outPath, '.'.join([baseNameSaveSummary, 'summaryData', 'tsv']))
    summaryData.to_csv(outName_summary, header = True, index = None, sep = '\t', float_format = '%.16f')
    print(outName_summary + ' saved')
    print('PROCESS FINISHED')


# Ask inputs to read the TSV file
## QUANTILES
dirPath_input = '/scratch/RDS-FMH-UNITI-RW/Q1'
inputName_input = 'cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.genes.samples.tsv'
# outPath_input = '/home/alba/phD/uniti_GBA_exomes.AB_GQ_DP/results/results/GBA_AF005/GBA_Q1'
outPath_input = '/scratch/RDS-FMH-UNITI-RW/Q1/'
inputAnalysesPath_input = '/scratch/RDS-FMH-UNITI-RW/Q1/type_analyses_AF005.tsv'
GBA_function(dirPath = dirPath_input, inputName = inputName_input, outPath = outPath_input, inputAnalysesPath = inputAnalysesPath_input)

print('ACABO')
