# Tiling deletion screen library design

A tiling deletion screen consists of a library where each mutation will delete a small portion of a protein of interest. The whole library can cover the full gene and help identify important domains for specific functions.

For example, in a screen consisiting of 20 residue deletions and a sliding window of 10 residues, the first mutation will delete resdiues 2 to 21, the second will delete residues 12 to 31, ...

![Tiling deletion drawing](https://github.com/aserracardona/CRISPR-scripts/blob/master/tiling_deletion_screen/tiling_del.png)

For each region to delete, we need a guide RNA targeting a PAM sequence close to the center of the region and a donor DNA containing both the upstream and downstream flanking sequences of the region to delete. 

## Library design
This library design is based on the CRISPR-Cas9 approach published in Guo et al. 2018 (https://doi.org/10.1038/nbt.4147). This library method assembles each guide and its correspondent donor sequences in the same plasmid. Each donor sequence consists of the 40bp upstream of the region to deleted + the 40bp downstream. 

The python script guide_donor_lib.py requires the gene name and accepts the optional arguments for the first and last amino acids to delete as well as how many residues to delete in each mutant, and the size of the sliding window between each mutant. It will download the gene sequence from the SGD database (www.yeastgenome.org) for the W303 background (for other backgrounds, edit line 45 of the code). 
The output will be an excel file containing the guide and donor sequences for all the different deletions in the selected gene. It also indicates the position of the NGG sequence closest to the center of the region to deleted. If there is no NGG sequence close enough (less than 11bp away by default) it won't retrive any guide sequence.

## Species
For now, it only works for *S. cerevisiae* genes. 

## Prerequisites
This script requieres [Biopython](https://biopython.org/) and [pandas](https://pandas.pydata.org/) installed. 

## Usage
python3 guide_donor_lib.py gene [-h] [-s START] [-e END] [-d DELETION] [-i INTERVAL]

positional arguments:
  gene         Gene name from S. cerevisiae

optional arguments:
  -h, --help   shows help message and exits
  -s START     Position of the first residue to delete in the tiling array. Second residue by default
  -e END       Position of the last residue to delete in the tiling array. Last amino acid of the protein by default
  -d DELETION  Number of residues to delete in each mutant. 20 by default
  -i INTERVAL  Interval of residues between on deletion and the next. 5 by default

Any feedback is welcome!
