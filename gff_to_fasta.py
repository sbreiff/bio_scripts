#!/usr/bin/env python

'''
Takes as input a gff file containing genes in a genome and the fasta file of that genome,
and outputs a fasta file of the genes specified in the gff file.

Example usage:
$ python gff_to_fasta.py genes.gff genome.fasta genes.fasta

'''

import csv
from Bio import SeqIO
from sys import argv
from SeqTools import rev_comp

gff_file = argv[1]
genome_file = argv[2]
out_file = argv[3]

# read in gff file, separate into rows
with open(gff_file, 'rb') as csvfile:
	reader = csv.reader(csvfile, delimiter = '\t')
	rows = list(reader)

# parse genome fasta, turn into list	
genome = list(SeqIO.parse(genome_file, "fasta"))

class Gene(object):

	def __init__(self, name, contig, strand):
		self.name = name.split(" ")[0]
		self.contig = contig
		self.coords = []
		self.seq = ""
		self.strand = strand
		
	def add_coords(self, coords):
		self.coords = coords
		
	def add_seq(self):
		for item in genome:
			if item.id == self.contig:
				seq = str(item.seq)
				break
		for item in self.coords:
			self.seq = self.seq + seq[int(item[0])-1:int(item[1])]
		if self.strand == "-":
			self.seq = rev_comp(self.seq)

# incorporate introns into sequence as lowercase
			
def create_fasta():
	with open(out_file, 'w') as outfile:
		gff_iterator = iter(rows)
		line = gff_iterator.next()
		while line:
			if line[2] == "gene":
				gene = Gene(line[8], line[0], line[6])
				coords = []
				line = gff_iterator.next()
				while line[2] != "gene":
					if line[2] == "CDS":
						coords.append([line[3], line[4]])
					try:
						line = gff_iterator.next()
					except StopIteration:
						break
				gene.add_coords(coords)
				gene.add_seq()
				outfile.write(">%s\n%s\n" % (gene.name, gene.seq))
			else:
				try:
					line = gff_iterator.next()
				except StopIteration:
					break
				
def test1(gff_list):
	gff_iterator = iter(gff_list)
	line = gff_iterator.next()
	while line:
		if line[2] == "gene":
			gene = Gene(line[8], line[0], line[6])
			coords = []
			line = gff_iterator.next()
			while line[2] != "gene":
				if line[2] == "CDS":
					coords.append([line[3], line[4]])
				try:
					line = gff_iterator.next()
				except StopIteration:
					break
			gene.add_coords(coords)
			gene.add_seq()
		else:
			try:
				line = gff_iterator.next()
			except StopIteration:
				break
	return gene

#gene1 = test1(rows)
#print gene1.name
#print gene1.contig
#print gene1.strand
#print gene1.coords
#print gene1.seq
			
create_fasta()
