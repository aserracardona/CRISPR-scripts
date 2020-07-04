# Tiling deletion screen library design

A tiling deletion screen consists of a library where each mutation will delete a small portion of a protein of interest. The whole library can cover the full gene and help identify important domains for specific functions.

For example, in a screen consisiting of 20 residue deletions and a sliding window of 10 residues, the first mutation will delete resdiues 2 to 21, the second will delete residues 12 to 31, ...

![Tiling deletion drawing](D:\Google Drive\Programing\Python\Projects\Guide_donor_library\guide_donor_seq)

For each region to delete, we need a guide RNA targeting the center of the region and a donor DNA containing both the upstream and downstream flanking sequences of the deleted region. 




- reference
- import biopython
- yeast but also others
