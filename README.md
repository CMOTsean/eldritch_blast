# REpeat COpy Number Estimation By Blast  (RECON-EBB)

This is a script I created to visualise and summarise the number of copies of a known repeat which were found in long reads sequenced by Nanopore MinION.

The utility of this script is to count the number of repeat copies in each read from a table of BLASTN results (using `-outfmt 6`), to identify exact number of copies and differences in copy number between alleles in cases of heterozygosity. However since long-reads are prone to error, some copies may have incorrect lengths in the BLAST output, or be missing altogether. A visualisation tool is included in the script to manually check the reads which are outliers in copy number.

### Requirements
Python 3

### Input
A tab separated file containing the results of a BLASTN search with:
 -  `query`:	the full sequence of the repeat +/- flanking DNA (recommended 1000bp either side). The flanking DNA is required for the script to identify reads containing a full, intact repeat region.
 -  `database`:	the long reads you are looking for the repeat in, converted from FASTQ to FASTA using `seqtk` or similar
 -  `outfmt`:	6  To make a tab separated output without headers.

### Output

####Summary

A tab separated file with 2 columns:
1. Column containing copy numbers up to the max number observed
2. Column containing the number of reads with each copy number

*Options*

`--read-names`
A file containing a list of each read with an intact repeat region and the corresponding read copy number.


**Reminder**

This script will only count the forward repeats in ODIRA cases. The total copy number of the amplified region, i.e Forward+Reverse, can be obtained by multiplying by 2 and subtracting 1 (2n-1).

####Plots

*Options*
`-p all` or `--plots all`
Generates a plot for each read in the input. Reads which have an intact, complete repeat region (i.e. the repeats are flanked by appropriate DNA) the plots are in black. Otherwise they are in red. 
Each hit in the read is plotted as a horizontal line using the co-ordinate of the hit in the read. For the ODIRA case, the hits are plotted on the top row for forward repeats, and bottom row for reverse repeats. For tandem cases, the reads are alternately placed in top or bottom rows to improve readability.
Numbers are placed above forward repeats to aid in counting from plot.

`-p flanked` or `--plots flanked`
Same as above but only the reads which have an intact, complete repeat region are plotted.

`-p None` or `--plots None` (DEFAULT)
No plots are generated.


###Usage tips

- This script is intended to estimate the copy number of tandem repeats; either simple tandem amplification or an ODIRA-like amplification.
- You must know the length of the repeat. For ODIRA cases, the forward length goes from leftmost inverted repeat to rightmost, inclusive. The reverse length is the part which is amplified in reverse, i.e. the sequence immediately flanked by the inner inverted repeats.
- You might notice some gaps in the plots generated from this script. This is due to errors in the reads preventing a successful BLAST hit at that specific repeat. This has no influence on the estimation of copy numbers as this uses the distance between the hits corresponding to the flanking DNA to estimate the number of repeats between them.
