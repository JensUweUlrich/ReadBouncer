from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import getopt,sys,re

prefix_length = 360

def usage():
    print("Usage: create360bpReads.py -i <input_fasta> -o <output_fasta> -p <prefix_length>")

try:
    options, remainder=getopt.getopt(sys.argv[1:], 'i:o:h:p')

except getopt.GetoptError as err:
    print(str(err))
    usage()
    sys.exit()

for opt, arg in options:
    if opt in ('-i'):
        input_file=arg
    if opt in ('-h'):
        usage()
        sys.exit()
    elif opt in ('-o'):
        output_file=arg

with open(output_file, 'w') as f_out:
    for record in SeqIO.parse(input_file, "fasta"):
        new_record = SeqRecord(seq=record.seq[:prefix_length], id=record.id)
        r=SeqIO.write(new_record, f_out, 'fasta')
        if r!=1: print('Error while writing sequence:  ' + record.id)