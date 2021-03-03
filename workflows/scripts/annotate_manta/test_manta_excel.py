import pandas as pd

df = pd.read_excel (r'/home/xshang/ws_manta-test/hanna_test_DNA70280_201023_AHCW35DSXY_somatic_mantaSV.vcf.xlsx')
#print (df[['ID','MATEID']].dropna())

# drops rows with NA in columns ID and/or MATEID
df = df.dropna(subset=['ID','MATEID'])
#print(df)

# this part of the script removes duplicates (translocations from both directions)

#print(df[['ID','MATEID']])
#row_indices_b = []
#row_indices_a = []
row_indices = []

for a in df['ID']:
    for b in df['MATEID']:
        if a == b:
            #print(a)
            #print(df[df['MATEID']==b].index.values)

            # get indices of duplicates
            row_index_b = df[df['MATEID']==b].index.values
            row_index_a = df[df['ID']==a].index.values
            #row_indices_b.append(row_index_b[0])
            #row_indices_a.append(row_index_a[0])

            # append indices as tuples for each matching pair
            row_indices.append((row_index_a[0], row_index_b[0]))
            #df.drop(row_index)
            #print(row_index[0])
            #df.drop(row_index[0], inplace = True)
#print(df.drop(row_index)[['ID','MATEID']])

#print(row_indices_b)
#print(row_indices_a)
#print(row_indices)

# prints pairs only once, with (a,b)
#print([tup for tup in row_indices if tup[0] < tup[1]])

# get pair only once (a,b) and remove the same pair (b,a)
a_b_pairs = [tup for tup in row_indices if tup[0] < tup[1]]

#print([i[1] for i in a_b_pairs])

# second element of each tuple = row indices to be removed
remove_indices = [i[1] for i in a_b_pairs]


df.drop(remove_indices, 0, inplace=True)
#print(df[['ID','MATEID']])

#print(df)


# this part of the script highlights genes from the gene list

# open the genelist
genelist = open("genelist.txt", "r")
#print(genelist.read())
genelist = genelist.readlines()

# make a list with genes from the genelist
gene_list = []
for gene in genelist:
    gene_list.append(gene.rstrip("\n"))
#print(gene_list)
pattern = '|'.join(gene_list)
# search for patterns of the genes in the genelist in certain column
print(df.loc[df['DEL/DUP Genecrossings'].str.contains(pattern, na=False)])




#print(df[['GeneInfo 1','GeneInfo 2', 'DEL/DUP Genecrossings']])


#print(df['GeneInfo 1'].str.contains("RRNAD"))

#print(genelist)
#df['Genelist'] = 'NO'
#for gene in genelist:
    #print(gene)
    #print(df['DEL/DUP Genecrossings'].str.contains(gene, na=False))
    #df.loc[df['DEL/DUP Genecrossings'].str.contains(gene, na=False), 'Genelist'] = "YES"

#print(df['Genelist'])



#df["Example"] = df["Name"].str.extract("(LB|RB)")[0] + " category"

#df['NEWcolumn']='NO'
#df.loc[df['COLUMN_to_Check'].str.contains(pattern), 'NEWcolumn'] = 'YES'


#Search_for_These_values = ['PEA', 'DEF', 'XYZ'] #creating list
#pattern = '|'.join(Search_for_These_values)