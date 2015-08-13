'''
Script for determining PCR product length

Usage: python PCRLength.py <fasta_file>

Script then supplies prompts:
Gene ID - supply name of sequence in fasta_file that primers are in
F primer - supply sequence of forward primer - spaces, uppercase, and lowercase allowed
R primer - supply sequence of reverse primer, 5' -> 3' - spaces, upper and lowercase allowed
Script then outputs expected PCR product length and starts prompt for next PCR product.

When finished, type "done" into Gene ID: prompt.

Future: perhaps create version where it reads in a tab-delimited file and outputs a list?
'''

from Bio import SeqIO
from sys import argv
from SeqTools import rev_comp
import readline

genes = list(SeqIO.parse(argv[1], "fasta"))

while True:
    gene_name = raw_input("\nGene ID: ")
    if gene_name == "done":
        break
    for gene in genes:
        if gene.id == gene_name:
            seq = str(gene.seq)
    f_primer = raw_input("Forward primer: ").upper().replace(' ', '')
    if f_primer not in seq:
        print "Forward primer not present in DNA sequence."
        continue
    r_primer = rev_comp(raw_input("Reverse primer: ").replace(' ', '')).upper()
    if r_primer not in seq:
        print "Reverse primer not present in DNA sequence."
        continue
    else:
        length = seq.index(r_primer) - seq.index(f_primer) + len(r_primer)
        # print seq.index(f_primer)
        # print seq.index(r_primer)
        print "PCR Product Length: %d bp" % length
