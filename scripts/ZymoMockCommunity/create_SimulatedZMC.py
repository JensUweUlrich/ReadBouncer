#!/usr/bin/python3

import shutil
import glob
import os
import subprocess

testdir = "testData/"
genomedir = "Genomes/"
simdir = testdir + "SimulatedReads/"
pbmodeldir = "/home/jens/software/pbsim2/data/"

if not os.path.exists(testdir):
	os.mkdir(testdir)
if not os.path.exists(simdir):
	os.mkdir(simdir)

for acc in ["80","85","90","95","98"]:
    accdir = simdir + "ACC" + acc + "HMM94/"
    if not os.path.exists(accdir):
	    os.mkdir(accdir)
    for filename in glob.glob(genomedir + "/*.fasta"):
        head, tail = os.path.split(filename)
        newname = os.path.splitext(tail)[0]
        if "Saccharomyces" in filename:
            pbsimcmd = ["pbsim", "--depth", "1", "--accuracy-mean", "0." + acc,"--length-min","360", "--length-max", "360", "--hmm_model", pbmodeldir + "R94.model", "--seed" , "20210617", filename]
        else:
            pbsimcmd = ["pbsim", "--depth", "20", "--accuracy-mean", "0." + acc,"--length-min","360", "--length-max", "360", "--hmm_model", pbmodeldir + "R94.model", "--seed" , "20210617", filename]
        subprocess.run(pbsimcmd)
        outfilename = accdir + newname + ".fastq"
        with open(outfilename, 'wb') as outfile:
            for filename in glob.glob("sd_*.fastq"):
                if filename == outfilename:
                    # don't want to copy the output into the output
                    continue
                with open(filename, 'rb') as readfile:
                    shutil.copyfileobj(readfile, outfile)
        
        for file in glob.glob("sd_*"):
        	os.remove(file)
