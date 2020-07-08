# Tiling deletion screen library design

A tiling deletion screen consists of a library where each mutation will delete a small portion of a protein of interest. The whole library can cover the full gene and help identify important domains for specific functions.

For example, in a screen consisiting of 20 residue deletions and a sliding window of 10 residues, the first mutation will delete resdiues 2 to 21, the second will delete residues 12 to 31, ...

![Tiling deletion drawing](https://github.com/aserracardona/CRISPR-scripts/blob/master/tiling_deletion_screen/tiling_del.png)

For each region to delete, we need a guide RNA targeting the center of the region and a donor DNA containing both the upstream and downstream flanking sequences of the deleted region. 

## Library design
This library design is based on the CRISPR-Cas9 approach published in Guo et al. 2018 (https://doi.org/10.1038/nbt.4147). This library method assembles each guide and its correspondent donor sequences in the same plasmid. Each donor sequence consists of the 40bp upstream of the region to deleted + the 40bp downstream. 

When executed, the python script guide_donor_lib.py will ask for the gene name and then it will download its sequence from the SGD database (www.yeastgenome.org). Next, it will ask which region of the gene we want to generate the different guide and donor sequences, how many residues to delete in each mutant, and the size of the sliding window between each mutant. 
The output will be an excel file containing the guide and donor sequences for all the different deletions in the selected gene. It also indicates the position of the NGG sequence closest to the center of the region to deleted. If there is no NGG sequence close enough (less than 11bp away by default) it won't retrive any guide sequence.

This script requieres [Biopython](https://biopython.org/) and [pandas](https://pandas.pydata.org/) installed. It only works for *S. cerevisiae* genes but I can adapt it to other organisms if somebody is interested. 

Any feedback is welcome!
