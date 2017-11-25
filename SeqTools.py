#!/usr/bin/env python

'''
This will be a set of functions that can be utilized generally for working with 
genetic data.

Currently contains:
rev_comp - find reverse complement of a DNA or RNA sequence.
genetic_code - returns dictionary containing the genetic code corresponding to
    the input number.
transcribe - replaces Ts with Us in a DNA sequence.
translate - translates a dna sequence into its (1-frame) protein sequence.
is_orf - checks whether a nucleotide sequence is an unbroken open reading frame.
count_nts - counts number of each nucleotide in a DNA sequence.
gc - calculates % GC content in a nucleotide sequence.
get_orfs - returns all unbroken open reading frames above a minimum length in a nucleotide sequence.
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
    rc_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'g':'c', 'c':'g', 't':'a', 'N':'N'}
    
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



def geneticCode(code):
    '''
    Given a numerical code, this function returns a dictionary 
    containing the appropriate genetic code.

    Genetic Codes:
    1 - Standard Code
    2 - Vertebrate Mitochondrial Code
    3 - Yeast Mitochondrial Code
    4 - Mold/Protozoan/Coelenterate Mitochondrial Code and 
    Mycoplasma/Spiroplasma Code
    5 - Invertebrate Mitochondrial Code
    6 - Ciliate/Dasycladacean/Hexamita Nuclear Code
    9 - Echinoderm/Flatworm Mitochondrial Code
    10 - Euplotid Nuclear Code
    11 - Bacterial/Archaeal/Plant Plastid Code
    12 - Alternative Yeast Nuclear Code
    13 - Ascidian Mitochondrial Code
    14 - Alternative Flatworm Mitochondrial Code
    15 - Blepharisma Nuclear Code
    16 - Chlorophycean Mitochondrial Code
    21 - Trematode Mitochondrial Code
    22 - Scenedesmus obliquus Mitochondrial Code
    23 - Thraustochytrium Mitochondrial Code
    24 - Pterobranchia Mitochondrial Code
    25 - Candidate Division SR1/Gracilibacteria Code   
    '''

    # standard genetic code
    genetic_code = {'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'TGC':'C',
            'TGT':'C', 'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'TTC':'F',
            'TTT':'F', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'CAC':'H',
            'CAT':'H', 'ATA':'I', 'ATC':'I', 'ATT':'I', 'AAA':'K', 'AAG':'K',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'TTA':'L', 'TTG':'L',
            'ATG':'M', 'AAC':'N', 'AAT':'N', 'CCA':'P', 'CCC':'P', 'CCG':'P', 
            'CCT':'P', 'CAA':'Q', 'CAG':'Q', 'AGA':'R', 'AGG':'R', 'CGA':'R',
            'CGC':'R', 'CGG':'R', 'CGT':'R', 'AGC':'S', 'AGT':'S', 'TCA':'S', 
            'TCC':'S', 'TCG':'S', 'TCT':'S', 'ACA':'T', 'ACC':'T', 'ACG':'T', 
            'ACT':'T', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'TGG':'W', 
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGA':'*'}
    
    # check for valid genetic code number
    if code not in range(1, 7) + range(9, 17) + range(21, 26):
        print "Not a valid genetic code."
        return None

    else: # alter genetic code dictionary depending on genetic code # provided
        if code in [2, 3, 4, 5, 9, 13, 14, 21, 24]:
            genetic_code['TGA'] = 'W'
        if code in [2, 3, 5, 13, 21]:
            genetic_code['ATA'] = 'M'
        if code in [5, 9, 14, 21, 24]:
            genetic_code['AGA'] = 'S'
        if code in [5, 9, 14, 21]:
            genetic_code['AGG'] = 'S'
        if code in [9, 14, 21]:
            genetic_code['AAA'] = 'N'
        if code in [16, 22]:
            genetic_code['TAG'] = 'L'
        if code == 2:
            genetic_code['AGA'] = '*'
            genetic_code['AGG'] = '*'
        elif code == 3:
            genetic_code['CTT'] = 'T'
            genetic_code['CTC'] = 'T'
            genetic_code['CTA'] = 'T'
            genetic_code['CTG'] = 'T'
            genetic_code['CGA'] = None
            genetic_code['CGC'] = None
        elif code == 6:
            genetic_code['TAA'] = 'Q'
            genetic_code['TAG'] = 'Q'
        elif code == 10: 
            genetic_code['TGA'] = 'C'
        elif code == 12:
            genetic_code['CTG'] = 'S'
        elif code == 13:
            genetic_code['AGG'] = 'G'
            genetic_code['AGA'] = 'G'
        elif code == 14:
            genetic_code['TAA'] = 'Y'
        elif code == 15:
            genetic_code['TAG'] = 'Q'
        elif code == 22:
            genetic_code['TCA'] = '*'
        elif code == 23: 
            genetic_code['TTA'] = '*'
        elif code == 24:
            genetic_code['AGG'] = 'K'
        elif code == 25:
            genetic_code['TGA'] = 'G'

    return genetic_code


def transcribe(dna):
    '''
    Transcribes a DNA sequence into its corresponding RNA sequence.
    Input is a DNA sequence which can be uppercase and/or lowercase.
    Output is an RNA sequence in the same case as the input.

    Example Usage:
    >>> transcribe('ATGTGA')
    AUGUGA
    >>> transcribe('ATGgtattaaatacacagTGA')
    AUGguauuaaauacacagUGA
    '''
    return dna.replace('T', 'U').replace('t', 'u')


def translate(dna, code = 1):
    '''
    Translates a DNA or RNA sequence to its corresponding protein sequence.
    Takes in a string representing a DNA/RNA sequence and an integer representing
    the genetic code/translation table to be used and returns a string 
    representing the protein sequence.
    
    Usage:   translate(dna, code = 1)
    
    Example Usage:
    >>> translate('ATGTGA')
    M*
    >>> translate('ATGTGA', 2)
    MW

    Parameters:
    dna - string comprising a dna or rna sequence - must only contain letters 
          A, C, G, and either T or U. Can be uppercase or lowercase.
    code - the number of the genetic code/translation table to be used. See
           below. If no code is given, standard code is used.
    '''

    # if dna parameter is an RNA sequence, convert to corresponding DNA sequence
    if 'U' in dna.upper() and 'T' not in dna.upper():
        dna = dna.upper().replace('U', 'T')

    geneticCodeDict = geneticCode(code)
    
    protein = ''
    for i in range(0, len(dna), 3): # translate each codon to its aa
        if 'N' in dna[i:i+3].upper():
            protein = protein + 'X'
        elif len(dna[i:]) >= 3:
            protein = protein + geneticCodeDict[dna[i:i+3].upper()]
    return protein

def is_orf(protein):
    if protein.startswith('M') and protein.endswith('*') and '*' not in protein[:-1]:
        return True
    return False


def count_nts(seq):
    '''
    '''
    if 'U' in seq and 'T' not in seq:
        return {'A': seq.upper().count('A'), 'C': seq.upper().count('C'), 
                'G': seq.upper().count('G'), 'T': seq.upper().count('T')}
    elif 'T' in seq and 'U' not in seq:
        return {'A': seq.upper().count('A'), 'C': seq.upper().count('C'), 
                'G': seq.upper().count('G'), 'U': seq.upper().count('U')}
    else:
        print "Sequence error: both Ts and Us in sequence"
        return


def gc(seq):
    '''
    '''
    counts = count_nts(seq)
    return counts['C'] + counts['G']/float(len(seq))

def get_orfs(seq, min_length = 100, genetic_code = 1, rc = False):
    '''
    '''
    #orfs = []
    #for i in range(len(seq) - min_length):
    #    if seq[i:i+3] == 'ATG':
    #        translation = translate(seq[i:], genetic_code)
    #        if '*' in translation:
    #            coding = seq[i:(3*translation.index('*'))+i+3]
    #            subseq = False
    #            if orfs:
    #                for orf in orfs:
    #                    if coding in orf:
    #                        subseq = True
    #            if not subseq:
    #                orfs.append(coding)
    coords = []
    for i in [0, 1, 2]:
        translation = translate(seq[i:], genetic_code)
        starts = [a for a, x in enumerate(translation) if x == 'M']
        coords = coords + [(start*3+i, translation.index('*', start)*3+i+3) for start in starts if '*' in translation[start:]]
    non_redundant = [(min([coord[0] for coord in [item for item in coords if item[1] == j]]), j) for j in list(set([x[1] for x in coords]))]
    non_redundant.sort()
    orfs = [seq[item[0]:item[1]] for item in non_redundant if item[1] - item[0] >= min_length]
    if rc:
        return orfs
    else:
        return orfs + get_orfs(rev_comp(seq), min_length, genetic_code, True)
