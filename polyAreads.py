'''
This is a script for analyzing 3'UTRs among polyA RNAseq reads. Inputs are a
fastq file of polyadenylated reads, a bed_file of these reads mapped to a set of
gene models, and a fasta file of the set of gene models they mapped to (not in 
that order). The output is a tab-delimited file which tells you the gene_id the 
read was mapped to, the number of bases between the stop codon and the polyA 
tail, and the sequence of the read with the portion mapped to the gene in 
uppercase and the rest lowercase.

This program assumes all reads are coding strand - ie, they all end in polyA,
and none begin with polyT. 

Usage:
$ python polyAreads.py <bed_file> <fastq_file> <genes_file> <output_file>

'''

from Bio import SeqIO; from sys import argv; import csv

with open(argv[1], 'rb') as bed_file:
    rows = list(csv.reader(bed_file, delimiter='\t'))

with open(argv[2], 'r') as fastq_file:
    reads = list(SeqIO.parse(fastq_file, 'fastq'))

with open(argv[3], 'r') as gene_file:
    genes = list(SeqIO.parse(gene_file, 'fasta'))

with open(argv[4], 'w') as out_file:
    for row in rows:
        for read in reads:
            # only reads that are same strand as gene
            if str(read.id) == row[3] and row[5] == "+":
                for gene in genes:
                    if row[0] == gene.id:
                        stop = len(str(gene.seq)) # 
                        break
                # make sure mapped portion of read extends to stop codon
                if int(row[2]) == stop:
                    idx_stop = int(row[2]) - int(row[1])
                    idx_tail = str(read.seq).index('AAAAAAAAAAAAAAAAAAAA')
                    out_file.write("%s\t%d\t%s%s\n" % (row[0], idx_tail - idx_stop, 
                                                       str(read.seq)[:idx_stop].upper(), 
                                                       str(read.seq)[idx_stop:].lower()))
                break
