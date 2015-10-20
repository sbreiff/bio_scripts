'''
Contains a couple versions of a fasta parser.

1. parse_fasta: A simple fasta parser

The parse_fasta function returns a list of Seq objects.
Each Seq object contains the sequence id, the sequence description,
and the sequence.

Example Usage (with example fasta file named 'samplefasta.fa'):
>>> from fasta_parser import parse_fasta
>>> records = parse_fasta('samplefasta.fa')
>>> for record in records:
...     print record.id
...
gene1
gene2
gene3
>>>

2. fasta_iterator: A fasta parser that creates a generator rather than loading
all the results into memory. Useful if fasta file is very large.

Used the same way as parse_fasta, except the contents can only be iterated
through once.

To be able to iterate through results more than once without creating the generator
again, you can convert the contents to a list (which will load all results into
memory at once):
>>> from fasta_parser import fasta_iterator
>>> records = list(fasta_iterator('samplefasta.fa'))

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
            
def fasta_iterator(file):
    fasta = open(file, 'r')
    line = fasta.readline()
    while True:
        if line.startswith(">"):
            text = line
            line = fasta.readline()
            while not line.startswith(">"):
                text += line
                line = fasta.readline()
                if not line:
                    break
            yield Seq(text)
        else:
            line = fasta.readline()
        if not line:
            return
