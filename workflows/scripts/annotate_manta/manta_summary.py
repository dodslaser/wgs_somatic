import pandas as pd
import re


def manta_summary(mantaSV_vcf, mantaSV_summary):
    df = pd.read_excel (mantaSV_vcf)

    # drops rows with NA in columns ID. MATEID can be NA if not a translocation
    df = df.dropna(subset=['ID'])


    # this part of the script removes duplicates (translocations from both directions)


    row_indices = []

    for a in df['ID']:
        for b in df['MATEID']:
            if a == b:

                row_index_b = df[df['MATEID']==b].index.values
                row_index_a = df[df['ID']==a].index.values


                # append indices as tuples for each matching pair
                row_indices.append((row_index_a[0], row_index_b[0]))



    # get pair only once (a,b) and remove the same pair (b,a)
    a_b_pairs = [tup for tup in row_indices if tup[0] < tup[1]]



    # second element of each tuple = row indices to be removed
    remove_indices = [i[1] for i in a_b_pairs]


    df.drop(remove_indices, 0, inplace=True)






    # this part of the script highlights genes from the gene list

    # open the genelist
    genelist = open("genelist.txt", "r")

    genelist = genelist.readlines()

    # make a list with genes from the genelist
    gene_list = []
    for gene in genelist:
        gene_list.append(gene.rstrip("\n"))
    # remove duplicates 
    gene_list = list(set(gene_list))

    # pattern to search for with an "or" between each gene
    pattern = '|'.join(gene_list)
    # search for patterns of the genes in the genelist in certain column
    column_patterns_genecrossings = df.loc[df['DEL/DUP Genecrossings'].str.contains(pattern, na=False)]['DEL/DUP Genecrossings']
    column_patterns_geneinfo1 = df.loc[df['GeneInfo 1'].str.contains(pattern, na=False)]['GeneInfo 1']
    column_patterns_geneinfo2 = df.loc[df['GeneInfo 2'].str.contains(pattern, na=False)]['GeneInfo 2']


    # create new column Genelist to contain genes from genelist that exist in certain columns
    df['Genelist'] = ""

    # function to find only whole words, otherwise it would find gene AR in XXARXX
    def find_only_whole_word(search_string, input_string):
        # Create a raw string with word boundaries from the user's input_string
        raw_search_string = r"\b" + search_string + r"\b"

        match_output = re.search(raw_search_string, input_string)

        no_match_was_found = ( match_output is None )
        if no_match_was_found:
            return False
        else:
            return True

    # function appending genes found to Genelist column
    def append_genes(col_patterns, col):
        for gene in gene_list:
            for row_of_genes in col_patterns:
                if gene in row_of_genes:
                    if find_only_whole_word(gene, row_of_genes) == True:
                        index_value = int(df[col==row_of_genes].index.values)
                        df.at[index_value, 'Genelist'] = df.at[index_value, 'Genelist'] + gene + ' '
        return df

    append_genes(column_patterns_genecrossings, df['DEL/DUP Genecrossings'])
    append_genes(column_patterns_geneinfo1, df['GeneInfo 1'])
    append_genes(column_patterns_geneinfo2, df['GeneInfo 2'])


    # remove duplicates in Genelist column
    for ind in df.index: 
        df['Genelist'][ind] = ' '.join(set(df['Genelist'][ind].split()))







    # this part of the script extracts allele frequencies from PR/SR


    # extract variants only supported by both PR and SR
    df = df.loc[df['FORMAT'] == 'PR:SR']



    # name of normal sample is in column name
    normal_sample = df.columns[31]

    # add columns for PR/SR for normal sample
    df[normal_sample + ':PR'] = ""
    df[normal_sample + ':PR-alt'] = ""
    df[normal_sample + ':SR'] = ""
    df[normal_sample + ':SR-alt'] = ""
    df['TOTAL alt (N)'] = ""
    df['TOTAL VAF (N)'] = ""

    # for row in column that contains PR/SR info for normal, add info to new columns
    for row in df.iloc[:, 31]:   
        row_index = int(df.loc[df[normal_sample] == row].index.values)
        PR = int(re.split('[,:]', row)[0])
        PR_alt = int(re.split('[,:]', row)[1])
        SR = int(re.split('[,:]', row)[2])
        SR_alt = int(re.split('[,:]', row)[3])
        df.at[row_index, normal_sample + ':PR'] = PR
        df.at[row_index, normal_sample + ':PR-alt'] = PR_alt
        df.at[row_index, normal_sample + ':SR'] = SR
        df.at[row_index, normal_sample + ':SR-alt'] =  SR_alt
        df.at[row_index, 'TOTAL alt (N)'] = PR_alt + SR_alt
        df.at[row_index, 'TOTAL VAF (N)'] = str(int(round(float(PR_alt + SR_alt) / (PR + PR_alt + SR + SR_alt) *100))) + '%'

    # name of tumor sample is in column name
    tumor_sample = df.columns[32]

    # add columns for PR/SR for tumor sample
    df[tumor_sample + ':PR'] = ""
    df[tumor_sample + ':PR-alt'] = ""
    df[tumor_sample + ':SR'] = ""
    df[tumor_sample + ':SR-alt'] = ""
    df['TOTAL alt (T)'] = ""
    df['TOTAL VAF (T)'] = ""

    # for row in column that contains PR/SR info for tumor, add info to new columns
    for row in df.iloc[:, 32]:
        row_index = int(df.loc[df[tumor_sample] == row].index.values)
        PR = int(re.split('[,:]', row)[0])
        PR_alt = int(re.split('[,:]', row)[1])
        SR = int(re.split('[,:]', row)[2])
        SR_alt = int(re.split('[,:]', row)[3])
        df.at[row_index, tumor_sample + ':PR'] = PR
        df.at[row_index, tumor_sample + ':PR-alt'] = PR_alt
        df.at[row_index, tumor_sample + ':SR'] = SR
        df.at[row_index, tumor_sample + ':SR-alt'] =  SR_alt
        df.at[row_index, 'TOTAL alt (T)'] = PR_alt + SR_alt
        df.at[row_index, 'TOTAL VAF (T)'] = str(int(round(float(PR_alt + SR_alt) / (PR + PR_alt + SR + SR_alt) *100))) + '%'



    df.to_excel(mantaSV_summary) 