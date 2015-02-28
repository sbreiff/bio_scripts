'''
Input: list of protein IDs
Output: fasta file containing all proteins specified in input file
Command line usage: $ python get_seqs.py <input_file>
'''

from Bio import SeqIO
from sys import argv

infile = argv[1]
outfile = infile[:-4] + "_seqs.fasta"

with open(infile, 'rb') as f:
	seqlist = f.readlines()
seqlist = [item.rstrip('\n') for item in seqlist]

proteome = list(SeqIO.parse("~/AllGeneModels_Translated.fasta",
		"fasta"))

with open(outfile, 'w') as out_f:
	for protein in proteome:
		if protein.id in seqlist:
                        out_f.write(">%s\n%s\n" % (protein.id, str(
                                protein.seq)))
