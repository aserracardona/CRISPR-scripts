from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import re
import urllib.request, json
import argparse

parser = argparse.ArgumentParser(description='This program asks for a budding yeast gene name, retrieves the sequence from the SGD database, '+ 
                                 'and generates the primers needed to design a deletion tiling array.')

parser.add_argument('gene', help='Gene name from S. cerevisiae')

parser.add_argument('-s', help='Position of the first residue to delete in the tiling array. Second residue by default', dest='start', type=int)

parser.add_argument('-e', help='Position of the last residue to delete in the tiling array. Last amino acid of the protein by default', dest='end', type=int)

parser.add_argument('-d', help='Number of residues to delete in each mutant. 20 by default', dest='deletion', type=int)

parser.add_argument('-i', help='Interval of residues between on deletion and the next. 5 by default', dest='interval', type=int)

args = parser.parse_args()

print(vars(args))

#variables from the arguments

gene_name = args.gene
        
starting_aa = args.start

ending_aa = args.end

del_aa = args.deletion

interval_aa = args.interval

max_cut_dist = 11 #This is the max distance from the center to deleted region to the cut site that we allow


def retrive_seq(name): #Function to download the gene sequence (+-1kb) from SGD
    with urllib.request.urlopen('https://www.yeastgenome.org/webservice/locus/PHO89/sequence_details') as url:
        data = json.loads(url.read().decode())
        
    for i in data['1kb']:
        if i['strain']['display_name'] == 'W303': #Select genetic background
            return i['residues']

sequence = retrive_seq(gene_name)

aa_num = int(((len(sequence)-2000)/3)-1) #Calculate aa length taking into account +1kb and -1kb

if starting_aa == None:
    starting_aa = 2
if ending_aa == None:
    ending_aa = aa_num
if del_aa == None:
    del_aa = 20
if interval_aa == None:
    interval_aa = 5

index_list = []

for i in range((starting_aa*3)+1000-3, (ending_aa*3)+1000, interval_aa*3):
    index_list.append(str(i+1) + '-' + str(i+del_aa*3))
    
gd_df = pd.DataFrame(columns=['aa_pos', 'ngg_pos', 'ccn_pos', 'best_cut',
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
                gd_df.loc[i, 'best_cut'] = 'No NGG sequence close to the center'
            else:
                guide = Seq(gd_df.loc[i, 'del_region'][gd_df.loc[i, 'ccn_pos']+3:gd_df.loc[i, 'ccn_pos']+24])
                gd_df.loc[i, 'best_cut'] = str(guide.reverse_complement())
        else:
            if abs(gd_df.loc[i, 'ngg_pos']-4-30) > max_cut_dist:
                gd_df.loc[i, 'best_cut'] = 'No NGG sequence close to the center'
            else:    
                gd_df.loc[i, 'best_cut'] = gd_df.loc[i, 'del_region'][gd_df.loc[i, 'ngg_pos']-22:gd_df.loc[i, 'ngg_pos']-1]
    except:
        pass
    
for i in gd_df.index:
    if gd_df.loc[i, 'ngg_pos'] == -1 and gd_df.loc[i, 'ccn_pos'] == -1:
        gd_df.loc[i, 'best_cut'] = 'No NGG sequence'
        
gd_df['donor'] = gd_df[['donor_fw', 'donor_rv']].apply(lambda x: ''.join(x), axis=1)

gd_df.to_excel('./' + gene_name + '_guide+donor.xlsx')

total = gd_df.count().loc['del_region']

if 'Too far' in gd_df.groupby('best_cut').count().index:
    too_far = gd_df.groupby('best_cut').count().loc['Too far', 'del_region']
else:
    too_far = 0
    
if 'No NGG sequence close to the center' in gd_df.groupby('best_cut').count().index:
    no_ngg = gd_df.groupby('best_cut').count().loc['No NGG sequence close to the center', 'del_region']
else:
    no_ngg = 0

print('Total number of deletions: ' + str(total) + '\nRegions without close NGG: ' +str(too_far) +'\nDeletions without NGG: ' +str(no_ngg))
