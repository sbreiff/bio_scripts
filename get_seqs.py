'''
Input: file containing list of protein IDs
Output: fasta file containing all proteins specified in input file
Command line usage: $ python get_seqs.py <input_file> <proteome_file>
'''

from Bio import SeqIO
from sys import argv

infile = argv[1]
outfile = infile[:-4] + "_seqs.fasta"

with open(infile, 'rb') as f:
	seqlist = [item.rstrip('\n') for item in f.readlines()]

with open(argv[2], 'r') as proteome_fasta:
        proteome = list(SeqIO.parse(proteome_fasta, "fasta"))

with open(outfile, 'w') as out_f:
	for protein in proteome:
		if protein.id in seqlist:
                        out_f.write(">%s\n%s\n" % (protein.id, str(
                                protein.seq)))
