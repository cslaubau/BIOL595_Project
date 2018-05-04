# BIOL595_Project

This project was submitted May 4th, 2018 by Curtis Slaubaugh, Rupesh Gaire, and Apostolia Topaloudi for BIOL 595

The script takes three inputs from the user: a gwas file, a gff file for the organism, and a fasta file of the
organism's cds. (Here the example inputs are 'Netblotch.txt', 'Hv_IBSC_PGSB_v2p37_genes.csv', and
'Hordeum_vulgare.Hv_IBSC_PGSB_v2.cds.all.fa')

A Manhattan plot can be generated from the GWAS file and -log(p) cutoff.
From these inputs, the user is asked for a -log(p) cutoff for which to isolate SNPs of interest from the GWAS file.
With the SNPs of interest, the program will generate a list of base pair ranges (with the range size determined by the
user, default 200k). The gff file will then be used to isolate coding named coding sequences for which these SNP regions
overlap. The named coding sequences are then retreived from the input fasta file and saved to a new key fasta file.

With the key fasta the user can perform a Blastx and InterPro search. The Blastx search will save a user defined number
of results for each cds. In this version of the project, the local search feature does not function. The fasta file is
split into smaller files containing a maximum of 10 coding sequences in order to reserve blast computing power.
The user can define an output file to which the blast results are saved.

The InterPro search takes the key fasta file, parses individual sequences, converts the sequence to amino acid
sequences, and searches the remote InterPro database. The returned xml data is then parsed and key information about
each result saved with respect to its sequence of origin. Sequences that do not begin with a start codon are excluded
from the searches.

It is not recommended that the Blast and Interpro search are performed in the same session (call of gui.py). They are
both remote processes that take a long time to process.

The code pieces used to generate the main program, gui.py, can be found in the folder 'Code bits'.