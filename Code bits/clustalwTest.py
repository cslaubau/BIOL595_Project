import os
from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo

clustalw_exe = r"./clustalw2"
clustalw_cline = ClustalwCommandline(clustalw_exe, infile="snpGenes.fasta")
assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
stdout, stderr = clustalw_cline()


tree = Phylo.read("snpGenes.dnd", "newick")

#Phylo.draw_ascii(tree)
tree.rooted = True
Phylo.draw(tree)