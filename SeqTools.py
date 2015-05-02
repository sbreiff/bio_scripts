#!/usr/bin/env python

'''
This will be a set of functions that can be utilized generally for working with 
genetic data.

Currently contains:
rev_comp - find reverse complement of a DNA or RNA sequence.

'''

def rev_comp(dna):
    '''
    Takes a string representing a DNA sequence as an argument, and returns the 
    reverse complement of that sequence. Input string can be uppercase or lowercase
    or a mix.
                                                                              
    Example usage:
    >>> rev_comp('tcAGCAgt')
    'acTGCTga'
    >>> rev_comp('agUCGAUcua')
    'uagAUCGAcu'

    '''

    # dictionary determining the reverse complement of each base
    rc_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'g':'c', 'c':'g', 't':'a'}
    
    # check whether sequence has Us and no Ts, in which case the dictionary is
    #     modified for RNA sequences   
    if 'U' in dna.upper() and 'T' not in dna.upper():
        rc_dict['A'] = 'U'
        rc_dict['U'] = 'A'
        rc_dict['a'] = 'u'
        rc_dict['u'] = 'a'

    for nt in dna:
        if nt not in rc_dict: # raise exception if invalid characters are used
            raise ValueError('Please provide a valid nucleotide sequence.')
            return None

    rc = [rc_dict[nt] for nt in dna]
    rc.reverse()
    return ''.join(rc)
