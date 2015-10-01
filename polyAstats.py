'''
This is a script for gathering statistics on 3'UTRs on a per gene basis. The
input is a tab-delimited file containing gene name, 3'UTR length, and read
sequence. The output is a tab-delimited file containing the gene name, the
mean 3'UTR length for that gene, the standard deviation, the minimum length,
the maximum length, the mode, and the number of data points.

Usage:
$ python polyAstats.py <input_file> <output_file>

***Actually may want to modify so it doesn't assume genes in input file are
sorted - make a dictionary instead***

'''

import csv; import numpy as np; from sys import argv

with open(argv[1], 'rb') as csvfile:
    rows = list(csv.reader(csvfile, delimiter='\t'))

# print min([int(row[1]) for row in rows])
# print len(set([row[0] for row in rows]))

lengths = []
gene = rows[0][0]
with open(argv[2], 'w') as outfile:
    outfile.write("gene\tmean length\tstdev\tmin\tmax\tmode\tn\n") # header
    for row in rows:
        if row[0] == gene:
            if int(row[1]) >= 0:
                lengths.append(int(row[1]))
            elif int(row[1]) in [-2, -1]:
            # reads with no UTR ending with TAA or TGA
                lengths.append(0)
        else:
            outfile.write("%s\t%g\t%s\t%d\t%d\t%d\t%d\n" % (
                    gene, np.mean(lengths), '{:.2f}'.format(np.std(lengths)),
                    min(lengths), max(lengths), np.bincount(lengths).argmax(),
                    len(lengths)))
            gene = row[0]
            lengths = [int(row[1]) if int(row[1]) >= 0 else 0]
    outfile.write("%s\t%g\t%s\t%d\t%d\t%d\n" % (gene, np.mean(lengths),
                                                '{:.2f}'.format(np.std(lengths)), 
                                                min(lengths), max(lengths), 
                                                np.bincount(lengths).argmax(), 
                                                len(lengths)))
