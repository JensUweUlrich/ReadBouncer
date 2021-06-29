from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import getopt,sys,re


def usage():
    print("Usage: splitByNs.py -i <input_scaffold_fasta> -o <output_contig_fasta>")

try:
    options, remainder=getopt.getopt(sys.argv[1:], 'i:o:h')

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
        sequence = ''.join(str(record.seq).strip())
        m=re.sub('[nN]+','\n',sequence).split('\n')
        for i in range(0,len(m)):
            if len(m[i]) > 100:
                new_record = SeqRecord(seq=Seq(m[i]), id=record.id + "_" + str(i))
                r=SeqIO.write(new_record, f_out, 'fasta')
                if r!=1: print('Error while writing sequence:  ' + record.id)



