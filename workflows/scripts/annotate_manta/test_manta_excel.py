import pandas as pd

df = pd.read_excel (r'/seqstore/webfolders/wgs/barncancer/hg38/DNA70280_test_manta/tumor/manta/DNA70280_201023_AHCW35DSXY_somatic_mantaSV.vcf.xlsx')
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
            row_index_b = df[df['MATEID']==b].index.values
            row_index_a = df[df['ID']==a].index.values
            #row_indices_b.append(row_index_b[0])
            #row_indices_a.append(row_index_a[0])
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

a_b_pairs = [tup for tup in row_indices if tup[0] < tup[1]]

# second element of each tuple = row indices to be removed
#print([i[1] for i in a_b_pairs])

remove_indices = [i[1] for i in a_b_pairs]


df.drop(remove_indices, 0, inplace=True)
print(df[['ID','MATEID']])



