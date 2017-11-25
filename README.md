### bio_scripts

Contains assorted scripts that I used to work with genetic and genomic data.


#### Contents:

1. SeqTools.py - collection of functions for dealing with sequence data.
    - rev_comp - reverse complements a nucleotide sequence.
    - genetic_code - returns dictionary containing the genetic code corresponding to the input number.
    - transcribe - replaces Ts in a DNA sequence with Us
    - translate - translates a DNA sequence into an amino acid sequence.
    - is_orf - checks whether a nucleotide sequence is an unbroken open reading frame.
    - count_nts - counts the number of each nucleotide in a sequence.
    - gc - calculates % GC content in a nucleotide sequence.
    - get_orfs - returns all the unbroken open reading frames above a minimum length within a nucleotide sequence.

2. gff_to_fasta.py - a script that takes in a genome fasta file and a gff file of genes, and outputs a fasta file of the genes in the gff.

3. get_seqs.py - a script that takes in a list of sequences and a fasta file that contains those sequences, and outputs a fasta file containing only the sequences in the list.

4. CodonUsage.py - computes codon usage in a sequence or a set of sequences.

5. PCRLength.py - a script that takes a fasta file as an argument, and prompts the user to supply a gene ID, a forward primer sequence, and a reverse primer sequence, and prints the length of the expected PCR product.

6. fasta_parser.py - for parsing fasta files

7. polyAreads.py - creates a table of information about polyadenylated RNA-Seq reads

8. polyAstats.py - takes table created from polyAreads.py and outputs a new table of 3'UTR stats for each gene in table.
