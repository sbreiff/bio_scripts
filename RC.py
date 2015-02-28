def reverse_complement(seq):
"""
Reverse complements a nucleotide sequence.
Input: nucleotide sequence (string) - ambiguous base codes not allowed
Output: reverse complement nucleotide sequence (string)
"""
	rc_dict = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A', 'a' : 't', 'c' : 'g',
		'g' : 'c', 't' : 'a'}
	rc = []

	try:
                for i in range(len(seq)):
			rc.append(rc_dict[seq[i]])
                rc.reverse()
                return ''.join(rc)

        # Raise a ValueError if a non-nucleotide character is present in string
        except KeyError:
                raise ValueError("Input is not a proper nucleotide sequence.")
