import subprocess
import shutil
import glob
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

testdir = "testData/"
tempdir = testdir + "temp/"
mockdir = testdir + "RealMock/"
genomedir = "Genomes/"
deepnanoreads = testdir + "allReads.DeepNanoW48.fasta"
fast5dir = "200416_MI2_Run20-087/fast5/"

if not os.path.exists(testdir):
	os.mkdir(testdir)
if not os.path.exists(tempdir):
	os.mkdir(tempdir)
if not os.path.exists(mockdir): 
   	os.mkdir(mockdir)

# start deepnano-blitz for basecalling
deepnanocmd = ["deepnano2_caller.py", "--output", deepnanoreads, "--directory", fast5dir, "--network-type", "48", "--output-format", "fasta", "--threads", "8"]
subprocess.run(deepnanocmd)

# creating reads to map with minimap2 => works
with open(testdir + "RealMockDataSet.fasta", 'w') as f_out:
    for record in SeqIO.parse(deepnanoreads, "fasta"):
        if len(record.seq) >= 2000:
            new_record = SeqRecord(seq=record.seq[360:], id=record.id)
            r=SeqIO.write(new_record, f_out, 'fasta')
            if r!=1: print('Error while writing sequence:  ' + record.id)

# make one reference file for mapping
outfilename = testdir + "AllZMCReferences.fasta"
with open(outfilename, 'wb') as outfile:
    for filename in glob.glob(genomedir + "*.fasta"):
        if filename == outfilename:
            # don't want to copy the output into the output
            continue
        with open(filename, 'rb') as readfile:
            shutil.copyfileobj(readfile, outfile)

# map reads against ZMC references
with open(os.devnull, 'w') as devnull:
    minimapcmd = ["minimap2","-ax","map-ont","-t","8",outfilename,testdir + "RealMockDataSet.fasta"]
    minimapoutput = subprocess.Popen(minimapcmd, stdout=subprocess.PIPE,stderr=devnull)
    samcmd = ["samtools","view", "-bh"]
    samoutput = subprocess.Popen(samcmd, stdin=minimapoutput.stdout, stdout=subprocess.PIPE, stderr=devnull)
    samsortcmd = ["samtools", "sort","-O","BAM"]
    c = subprocess.Popen(samsortcmd ,stdin=samoutput.stdout, stdout=open(tempdir + "RealZMC.bam", 'w'), stderr=devnull)
    minimapoutput.stdout.close()
    samoutput.stdout.close()
    c.communicate()
    subprocess.run(["samtools", "index",tempdir + "RealZMC.bam"], stderr=devnull)
    
    for filename in glob.glob(genomedir + "/*.fasta"):
        head, tail = os.path.split(filename)
        newname = os.path.splitext(tail)[0]
        samviewcmd = ["samtools","view","-q","30","-bh",tempdir + "RealZMC.bam"]
        for record in SeqIO.parse(filename, "fasta"):
            samviewcmd.append(record.id)

        samviewout = subprocess.Popen(samviewcmd, stdout=subprocess.PIPE, stderr=devnull)
        proc = subprocess.Popen(samsortcmd ,stdin=samviewout.stdout, stdout=open(tempdir + newname + ".bam", 'w'), stderr=devnull)
        samviewout.stdout.close()
        proc.communicate()
        samfastacmd = ["samtools","fasta",tempdir + newname + ".bam"]
        fproc = subprocess.Popen(samfastacmd , stdout=open(tempdir + newname + ".fasta", 'w'), stderr=devnull)
        fproc.communicate()
        readids = []
        for record in SeqIO.parse(tempdir + newname + ".fasta", "fasta"):
            readids.append(record.id)
        
        with open(mockdir + newname + ".fasta", 'w') as f_out:
            for record in SeqIO.parse(deepnanoreads, "fasta"):
                if record.id in readids:
                    new_record = SeqRecord(seq=record.seq[:360], id=record.id)
                    r=SeqIO.write(new_record, f_out, 'fasta')
                    if r!=1: print('Error while writing sequence:  ' + record.id)

