import subprocess
import shutil
import glob
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# start deepnano-blitz for basecalling
#subprocess.run(["ls", "-l"])

# creating reads to map with minimap2
with open("../testData/RealMockDataSet.fasta", 'w') as f_out:
    for record in SeqIO.parse("../testData/allReads.DeepNanoW48.fasta", "fasta"):
        if len(record.seq) >= 2000:
            new_record = SeqRecord(seq=record.seq[720:], id=record.id)
            r=SeqIO.write(new_record, f_out, 'fasta')
            if r!=1: print('Error while writing sequence:  ' + record.id)

# make one reference file for mapping
outfilename = "../testData/Genomes/AllZMCReferences.fasta"
with open(outfilename, 'wb') as outfile:
    for filename in glob.glob('../testData/Genomes/*.fasta'):
        if filename == outfilename:
            # don't want to copy the output into the output
            continue
        with open(filename, 'rb') as readfile:
            shutil.copyfileobj(readfile, outfile)

tempdir = "../testData/temp/"
os.mkdir(tempdir)
# map reads against ZMC references
with open(os.devnull, 'w') as devnull:
#    minimapcmd = ["minimap2","-ax","map-ont","-t","8",outfilename,"../testData/RealMockDataSet.fasta"]
#    minimapoutput = subprocess.Popen(minimapcmd, stdout=subprocess.PIPE,stderr=devnull)
#    samcmd = ["samtools","view", "-bh"]
#    samoutput = subprocess.Popen(samcmd, stdin=minimapoutput.stdout, stdout=subprocess.PIPE, stderr=devnull)
#    samsortcmd = ["samtools", "sort","-O","BAM"]
#    c = subprocess.Popen(samsortcmd ,stdin=samoutput, stdout=open(tempdir + "RealZMC.bam", 'w'), stderr=devnull)
#    minimapoutput.stdout.close()
#    samoutput.stdout.close()
#    c.communicate()
    os.remove(outfilename)
    mockdir ="../testData/RealMock/" 
    os.mkdir(mockdir)
    for filename in glob.glob('../testData/Genomes/*.fasta'):
        head, tail = os.path.split(filename)
        newname = os.path.splitext(tail)[0]
        refs = []
        for record in SeqIO.parse(filename, "fasta"):
            refs.append(record.id)
        samviewcmd = ["samtools","view","-bh","-q","30",tempdir + "RealZMC.bam"]
        samviewcmd.append(refs)
        samviewout = subprocess.Popen(samviewcmd, stdout=subprocess.PIPE, stderr=devnull)
        proc = subprocess.Popen(samsortcmd ,stdin=samviewout, stdout=open(tempdir + newname + ".bam", 'w'), stderr=devnull)
        samviewout.stdout.close()
        proc.communicate()
        samfastacmd = ["samtools","fasta",tempdir + newname + ".bam"]
        fproc = subprocess.Popen(samfastacmd , stdout=open(tempdir + newname + ".fasta", 'w'), stderr=devnull)
        fproc.communicate()
        readids = []
        for record in SeqIO.parse(tempdir + newname + ".fasta", "fasta"):
            readids.append(record.id)
        
        with open(mockdir + newname + ".fasta", 'w') as f_out:
            for record in SeqIO.parse("../testData/allReads.DeepNanoW48.fasta", "fasta"):
                if record.id in readids:
                    r=SeqIO.write(record, f_out, 'fasta')
                    if r!=1: print('Error while writing sequence:  ' + record.id)

