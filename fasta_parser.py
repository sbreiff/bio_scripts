'''
A simple fasta parser

The parse_fasta function returns a list of Seq objects.
Each Seq object contains the sequence id, the sequence description,
and the sequence.

Example Usage (with example fasta file named 'samplefasta.fa'):
>>> from fasta_parser import parse_fasta
>>> records = parse_fasta('samplefasta.fa')

For info on fasta format, see https://en.wikipedia.org/wiki/FASTA_format#Format

'''

class Seq:
    def __init__(self, text):
        self.id = text.split()[0][1:]
        self.description = text.split('\n')[0][1:]
        self.seq = ''.join(text.split('\n')[1:])

def parse_fasta(file):
    with open(file, 'r') as fasta_file:
        lines = fasta_file.readlines()
    record = ''
    records = []
    for line in lines:
        if line.startswith('>') and not record:
            record = line
        elif not line.startswith('>') and record:
            record += line
        else:
            records.append(Seq(record))
            record = line
    records.append(Seq(record))
    return records
            
    
