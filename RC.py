"""
Reverse complements a nucleotide sequence.
Input: nucleotide sequence (string)
Output: reverse complement nucleotide sequence (string)
"""

def reverse_complement(seq):
	dict = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A', 'a' : 't', 'c' : 'g',
		'g' : 'c', 't' : 'a'}
	rc = []
	for i in range(len(seq)):
		if seq[i] in dict:
			rc.append(dict[seq[i]])
		else: 
			print "Input is not a proper nucleotide sequence."
			return "Input is not a proper nucleotide sequence."
	rc.reverse()
	with open("rc_out.txt", 'w') as f:
		f.write(''.join(rc))
	return ''.join(rc)