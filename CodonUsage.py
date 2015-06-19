from SeqTools import geneticCode
from itertools import product
from sys import argv

def codonUsage(fasta_file, out_file, code = 1):
    '''
    Takes in a fasta file containing a set of CDSs, and outputs a 
    tab-delimited file describing the codon usage of that set of genes into
    the specified output file. Uses the standard genetic code.
    '''

    # create dictionary to store how often each codon appears
    codon_dict = {} 
    for codon in [''.join(item) for item in list(product('ACGT', repeat=3))]:
        codon_dict[codon] = 0
    
    # dictionary that stores which amino acid each codon codes for
    gen_code = geneticCode(code)

    genes = ""
    with open(fasta_file, 'r') as fastafile: 
        # concatenate all gene seqs into 1 string
        for line in fastafile.readlines():
            if not line.startswith(">"):
                genes += line.rstrip('\n')
        
    for i in range(len(genes)/3): # count codons across entire combined sequence
        codon = genes[3*i:3*(i+1)]
        if 'N' in codon:
            continue
        codon_dict[codon] += 1

    grandtotal = float(sum(codon_dict.values())) # computes total codons in set
    
    with open(out_file, 'w') as outfile:
        outfile.write("AA\tcodon\t#\t% aa\t% all\n")
        outfile.write("--\t-----\t-\t----\t-----\n\n")
        for aa in 'ACDEFGHIKLMNPQRSTVWY*':
            # calculates total codons for each amino acid
            total = float(sum([codon_dict[codon] for codon in gen_code if 
                    gen_code[codon] == aa]))
            # create list of codons for that amino acid
            for codon in [cdn for cdn in gen_code if gen_code[cdn] == aa]:
                if total == 0: # avoid division by zero
                    t_pct = "0.00%"
                    gt_pct = "0.00%"
                else:
                    t_pct = "{0:.2f}%".format(codon_dict[codon]/total*100)
                    gt_pct = "{0:.2f}%".format(codon_dict[codon]/grandtotal*100)
                # for each codon, output its amino acid, the codon, no. of
                # appearances, % for that amino acid, and % of all codons
                outfile.write("%s\t%s\t%d\t%s\t%s\n" % (aa, codon, 
                            codon_dict[codon], t_pct, gt_pct))
            outfile.write('\n') # add black line after each amino acid


codonUsage(argv[1], argv[2]) # call function
