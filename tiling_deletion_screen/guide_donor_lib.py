from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import re
import urllib.request, json

gene_name = input('Type a gene name: ') #Retrieve the gene name from the window

def retrive_seq(name): #Function to download the gene sequence (+-1kb) from SGD
    with urllib.request.urlopen('https://www.yeastgenome.org/webservice/locus/' + gene_name + '/sequence_details') as url:
        data = json.loads(url.read().decode())
        
    for i in data['1kb']:
        if i['strain']['display_name'] == 'W303': #Select genetic background
            return i['residues']
        
sequence = retrive_seq(gene_name)

aa_num = int(((len(sequence)-2000)/3)-1) #Calculate aa length taking into account +1kb and -1kb

starting_aa = int(input('First residue to delete (2 - '+str(aa_num)+'): '))

ending_aa = int(input('Last residue to delete: (2 - '+str(aa_num)+'): '))

del_aa = int(input('How many residues to delete in each mutant: '))

interval_aa = int(input('Sliding window of residues between one deletion and the next: '))

max_cut_dist = 11 #This is the max distance from the center to deleted region to the cut site that we allow

index_list = []

for i in range((starting_aa*3)+1000-3, (ending_aa*3)+1000, interval_aa*3):
    index_list.append(str(i+1) + '-' + str(i+del_aa*3))
    
gd_df = pd.DataFrame(columns=['aa_pos', 'ngg_pos', 'ccn_pos', 'best_sgRNA',
                              'del_region', 'donor_fw', 'donor_rv'], index=index_list)

#determine the region to delete and the 40bp upstream donor and the 40bp downstream donor
for i in range((starting_aa*3)+1000-2, (ending_aa*3)+1000, interval_aa*3):
    gd_df.loc[str(i) + '-' + str(i-1+del_aa*3), 'aa_pos'] = str((i-1000+2)/3) + '-' + str((i-1000-1)/3+del_aa)
    gd_df.loc[str(i) + '-' + str(i-1+del_aa*3), 'del_region'] = sequence[i-1:i-1+del_aa*3]
    gd_df.loc[str(i) + '-' + str(i-1+del_aa*3), 'donor_fw'] = sequence[i-41:i-1]
    gd_df.loc[str(i) + '-' + str(i-1+del_aa*3), 'donor_rv'] = sequence[i-1+del_aa*3:i+del_aa*3+39]
    
#identify NGG and CCN sequences and select the one with the shortest cutting distance from the center of the deleted region
for i in gd_df.index:
    try:
        gd_df.loc[i, 'ngg_pos'] = min([m.start()-4-30 for m in re.finditer('GG', gd_df.loc[i,'del_region'])], key=abs) + 34
    except:
        pass

    try:
        gd_df.loc[i, 'ccn_pos'] = min([m.start()+6-30 for m in re.finditer('CC', gd_df.loc[i,'del_region'])], key=abs) + 24
    except:
        pass  
    
#select the closest cutting site to the center (NGG or CCN) and store the guide RNA for that site
gd_df = gd_df.fillna(-1) #otherwise the NaN values in ngg_pos and ccn_pos give error

for i in gd_df.index:
    try:
        if abs(gd_df.loc[i, 'ngg_pos']-4-30) > abs(gd_df.loc[i, 'ccn_pos']+6-30):
            if abs(gd_df.loc[i, 'ccn_pos']+6-30) > max_cut_dist:
                gd_df.loc[i, 'best_sgRNA'] = 'No NGG sequence close to the center'
            else:
                guide = Seq(gd_df.loc[i, 'del_region'][gd_df.loc[i, 'ccn_pos']+3:gd_df.loc[i, 'ccn_pos']+24])
                gd_df.loc[i, 'best_sgRNA'] = str(guide.reverse_complement())
        else:
            if abs(gd_df.loc[i, 'ngg_pos']-4-30) > max_cut_dist:
                gd_df.loc[i, 'best_sgRNA'] = 'No NGG sequence close to the center'
            else:    
                gd_df.loc[i, 'best_sgRNA'] = gd_df.loc[i, 'del_region'][gd_df.loc[i, 'ngg_pos']-22:gd_df.loc[i, 'ngg_pos']-1]
    except:
        pass
    
for i in gd_df.index:
    if gd_df.loc[i, 'ngg_pos'] == -1 and gd_df.loc[i, 'ccn_pos'] == -1:
        gd_df.loc[i, 'best_sgRNA'] = 'No NGG sequence'
        
gd_df['donor'] = gd_df[['donor_fw', 'donor_rv']].apply(lambda x: ''.join(x), axis=1)

gd_df.to_excel('./' + gene_name + '_guide+donor.xlsx')

total = gd_df.count().loc['del_region']

if 'Too far' in gd_df.groupby('best_sgRNA').count().index:
    too_far = gd_df.groupby('best_sgRNA').count().loc['Too far', 'del_region']
else:
    too_far = 0

if 'No NGG sequence close to the center' in gd_df.groupby('best_sgRNA').count().index:
    no_ngg = gd_df.groupby('best_sgRNA').count().loc['No NGG sequence close to the center', 'del_region']
else:
    no_ngg = 0


print('\n### '+gene_name.upper()+' Tiling Deletion Library ###\n'+
    'Total number of deletions: ' + str(total) + '\nRegions without close NGG: ' +str(too_far) +'\nRegions without NGG: ' +str(no_ngg)+'\n'+
    'File saved as ' + gene_name.upper() + '_guide+donor.xlsx')