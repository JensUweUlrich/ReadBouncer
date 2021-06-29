#!/usr/bin/python3

import shutil
import glob
import os
import subprocess
import time
import math
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import mappy


testdir = "testData/"
genomedir = "Genomes/"
simdir = testdir + "SimulatedReads/"
zmc360 = testdir + "RealMock_360bp/"
zmc180 = testdir + "RealMock_180bp/"
mockdir = testdir + "RealMock/"



monipath = "/home/jens/tools/spumoni/build/moni"
analyzepath = "/home/jens/tools/spumoni/analysis/analyze_pml.py"

resultdir = "results/"
spumonidir = resultdir + "spumoni/"
indexdir = spumonidir + "index/"
outdir = spumonidir + "output/"
realpmldir = outdir + "RealZMC_360bp_PML/"


minimapdir = resultdir + "minimap2/"
minioutdir = minimapdir + "output/"
realminidir = minioutdir + "RealZMC_180bp/"


if not os.path.exists(resultdir):
	os.mkdir(resultdir)
if not os.path.exists(spumonidir):
	os.mkdir(spumonidir)
if not os.path.exists(indexdir):
	os.mkdir(indexdir)
if not os.path.exists(outdir):
	os.mkdir(outdir)
if not os.path.exists(realpmldir):
	os.mkdir(realpmldir)
if not os.path.exists(zmc360):
	os.mkdir(zmc360)
if not os.path.exists(minimapdir):
	os.mkdir(minimapdir)
if not os.path.exists(minioutdir):
	os.mkdir(minioutdir)
if not os.path.exists(realminidir):
	os.mkdir(realminidir)	


def create_ZMC360bp():
    for filename in glob.glob(mockdir + "*.fasta"):
        head, tail = os.path.split(filename)
        with open(zmc360 + tail, 'w') as f_out:
            for record in SeqIO.parse(filename, "fasta"):
                new_record = SeqRecord(seq=record.seq[:360], id=record.id)
                r=SeqIO.write(new_record, f_out, 'fasta')
                if r!=1: print('Error while writing sequence:  ' + record.id)

# create bacterial references with their reverse complement sequences and yeast-only reference file
# add forward and reverse complement of every reference sequence to the files
# spumoni needs both for index building
def prepare_seqs_for_indexing():

    # make one reference file for mapping
    outfilename = testdir + "ZMCBacReferences.fasta"
    with open(outfilename, 'w') as f_out:
        for filename in glob.glob(genomedir + "*.fasta"):
            head, tail = os.path.split(filename)
            if tail.startswith("Saccharomyces"):
                continue
            for record in SeqIO.parse(filename, "fasta"):
                #new_record = SeqRecord(seq=record.seq[0:], id=record.id)
                r=SeqIO.write(record, f_out, 'fasta')
                if r!=1: print('Error while writing sequence:  ' + record.id)
                new_record = SeqRecord(seq=record.seq.reverse_complement(), id=record.id + "_revComp")
                r=SeqIO.write(new_record, f_out, 'fasta')
                if r!=1: print('Error while writing sequence:  ' + new_record.id)
     

def create_spumoni_indexes():
    #build positive index for bacterial references
    monibuildcmd = ["python3",monipath, "build", "-r", testdir + "ZMCBacReferences.fasta", "-f", "--spumoni"]
    subprocess.run(monibuildcmd)
    for file in glob.glob(testdir + "ZMCBacReferences.fasta.*"):
        shutil.move(file, indexdir)

    # build null index for bacterial references
    with open(testdir + "ZMCBacReferences_null.fasta", 'w') as f_out:
        for record in SeqIO.parse(testdir + "ZMCBacReferences.fasta", "fasta"):
            new_record = SeqRecord(seq=record.seq[::-1], id=record.id)
            r=SeqIO.write(new_record, f_out, 'fasta')
            if r!=1: print('Error while writing sequence:  ' + record.id)

    monibuildcmd = ["python3",monipath, "build", "-r", testdir + "ZMCBacReferences_null.fasta", "-f", "--spumoni"]
    subprocess.run(monibuildcmd)
    for file in glob.glob(testdir + "ZMCBacReferences_null.fasta.*"):
        shutil.move(file, indexdir)
        

def compute_spumoni_pmls_realZMC():

    tp = 0
    fp = 0
    tn = 0
    fn = 0
    complete_elapsed = 0.0
    org_counter = 0
    nr_allreads = 0
    with open(spumonidir +  "pml_analysis_ZMCReal_360bp.txt", 'w') as resultfile:
        for filename in glob.glob(zmc360 + "*.fasta"):
            head, tail = os.path.split(filename)
            newname = os.path.splitext(tail)[0]
            org_counter += 1
            start = time.perf_counter()
            # pseudo-ms for positive index of bacterial references
            pseudomscmd = ["python3", monipath, "pseudo-ms", "-i", indexdir + "ZMCBacReferences.fasta", "-p", filename]
            subprocess.run(pseudomscmd)
            stop = time.perf_counter()
            complete_elapsed += stop - start
            # pseudo-ms for null index of bacterial references
            pseudomscmd = ["python3", monipath, "pseudo-ms", "-i", indexdir + "ZMCBacReferences_null.fasta", "-p", filename]
            subprocess.run(pseudomscmd)
            for file in glob.glob(filename + "_*"):
                shutil.move(file, realpmldir)
            os.remove(filename + ".pseudo_ms.log")           
                        
            analyzecmd = ["python3", analyzepath, "-p", realpmldir + tail + "_ZMCBacReferences.fasta.pseudo_lengths", "-n", realpmldir + tail + "_ZMCBacReferences_null.fasta.pseudo_lengths", "-r", "90"]
            proc = subprocess.Popen(analyzecmd , stdout=open(realpmldir + newname + "_pml_analysis.txt", 'w'))
            proc.communicate()    
            
            nrmatches = 0
            nrreads = 0
            with open(realpmldir + newname + "_pml_analysis.txt", 'r') as f:
                readid = ""
                found = False
                for line in f:
                    match = line.split()
                    if int(match[2]) > 360:
                        continue
                    if readid != match[0]:
                        readid = match[0]
                        found = False
                        nrreads += 1
                    if match[4] == "found" and found == False:
                        found = True
                        nrmatches += 1
                        
            resultfile.write(newname + "\tZMCBacReferences\t" + str(nrmatches) + "\t" + str(nrreads) + "\n")
            nr_allreads += nrreads
            if not newname.startswith("Saccharomyces"):
                tp += int(nrmatches)
                fn += int(nrreads) - int(nrmatches)
            else:
                tn = int(nrreads) - int(nrmatches)
                fp = int(nrmatches)
            
            
            
        accuracy = float(tp + tn)/float(tp + tn + fp + fn)
        precision = float(tp)/float(tp + fp)
        recall = float(tp)/float(tp + fn)
        specificity = float(tn)/float(tn + fp)
        f1score = float(2*tp)/float(2*tp + fp + fn)
        mcc = float(tp*tn - fp*fn)/math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
        
        resultfile.write("Accuracy\t: " + str(accuracy) + "\n")
        resultfile.write("Precision\t: " + str(precision) + "\n")
        resultfile.write("Recall\t\t: " + str(recall) + "\n")
        resultfile.write("Specificity\t: " + str(specificity) + "\n")
        resultfile.write("F1-Score\t: " + str(f1score) + "\n")
        resultfile.write("MCC\t: " + str(mcc) + "\n")
        resultfile.write("Throughput (bp/sec)\t: " + str(float(nr_allreads*360)/complete_elapsed) + "\n")
    
def compute_spumoni_pmls_simulated():

    for sdir in glob.glob(simdir + "ACC*"):
        simdirpath, simdirname = os.path.split(sdir)
        accoutputdir = outdir + simdirname + "_360bp/"
        if not os.path.exists(accoutputdir):
	        os.mkdir(accoutputdir)
	        
        tp = 0
        fp = 0
        tn = 0
        fn = 0
  
        with open(spumonidir +  "pml_analysis_" + simdirname +"_360bp.txt", 'w') as resultfile:
            for filename in glob.glob(sdir + "/*.fastq"):
                head, tail = os.path.split(filename)
                newname = os.path.splitext(tail)[0]

                # pseudo-ms for positive index of bacterial references
                pseudomscmd = ["python3", monipath, "pseudo-ms", "-i", indexdir + "ZMCBacReferences.fasta", "-p", filename, "-t", "8"]
                subprocess.run(pseudomscmd)
                # pseudo-ms for null index of bacterial references
                pseudomscmd = ["python3", monipath, "pseudo-ms", "-i", indexdir + "ZMCBacReferences_null.fasta", "-p", filename, "-t", "8"]
                subprocess.run(pseudomscmd)
                for file in glob.glob(filename + "_*"):
                    shutil.move(file, accoutputdir)
                os.remove(filename + ".pseudo_ms.log")           
                        
                analyzecmd = ["python3", analyzepath, "-p", accoutputdir + tail + "_ZMCBacReferences.fasta.pseudo_lengths", "-n", accoutputdir + tail + "_ZMCBacReferences_null.fasta.pseudo_lengths", "-r", "90"]
                proc = subprocess.Popen(analyzecmd , stdout=open(accoutputdir + newname + "_pml_analysis.txt", 'w'))
                proc.communicate()
                nrmatches = 0
                nrreads = 0
                with open(accoutputdir + newname + "_pml_analysis.txt", 'r') as f:
                    readid = ""
                    found = False
                    for line in f:
                        match = line.split()
                        if int(match[2]) > 360:
                            continue
                        if readid != match[0]:
                            readid = match[0]
                            found = False
                            nrreads += 1
                        if match[4] == "found" and found == False:
                            found = True
                            nrmatches += 1
                        
                resultfile.write(newname + "\tZMCBacReferences\t" + str(nrmatches) + "\t" + str(nrreads) + "\n")
                if not newname.startswith("Saccharomyces"):
                    tp += int(nrmatches)
                    fn += int(nrreads) - int(nrmatches)
                else:
                    tn = int(nrreads) - int(nrmatches)
                    fp = int(nrmatches)
            
                
            
            accuracy = float(tp + tn)/float(tp + tn + fp + fn)
            precision = float(tp)/float(tp + fp)
            recall = float(tp)/float(tp + fn)
            specificity = float(tn)/float(tn + fp)
            f1score = float(2*tp)/float(2*tp + fp + fn)
            mcc = float(tp*tn - fp*fn)/math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
        
            resultfile.write("Accuracy\t: " + str(accuracy) + "\n")
            resultfile.write("Precision\t: " + str(precision) + "\n")
            resultfile.write("Recall\t\t: " + str(recall) + "\n")
            resultfile.write("Specificity\t: " + str(specificity) + "\n")
            resultfile.write("F1-Score\t: " + str(f1score) + "\n")
            resultfile.write("MCC\t: " + str(mcc) + "\n")

def compute_minimap_matches_realZMC():

    # make one reference file for mapping
    bacreffile = testdir + "ZMCForwardBacReferences.fasta"
    with open(bacreffile, 'wb') as outfile:
        for filename in glob.glob(genomedir + "*.fasta"):
            head, tail = os.path.split(filename)
            if tail.startswith("Saccharomyces"):
                continue
            with open(filename, 'rb') as readfile:
                shutil.copyfileobj(readfile, outfile)

    mapper = mappy.Aligner(bacreffile, preset="map-ont",n_threads=1)

    tp = 0
    fp = 0
    tn = 0
    fn = 0
    complete_elapsed = 0.0
    with open(minimapdir +  "minimap2_matching_analysis_ZMCReal_360bp_2.txt", 'w') as resultfile:
        start = time.perf_counter()
        for filename in glob.glob(zmc360 + "*.fasta"):
            head, tail = os.path.split(filename)
            newname = os.path.splitext(tail)[0]
            
            # map reads
            mapped = 0
            nrreads = 0
        
            for record in SeqIO.parse(filename, "fasta"):
                nrreads += 1
                for hit in mapper.map(record.seq):
                    mapped += 1
                    break
        
            
            nr_allreads = nrreads         
            resultfile.write(newname + "\tZMCForwardBacReferences\t" + str(mapped) + "\t" + str(nrreads) + "\n")
            if not newname.startswith("Saccharomyces"):
                tp += int(mapped)
                fn += int(nrreads) - int(mapped)
            else:
                tn = int(nrreads) - int(mapped)
                fp = int(mapped)
        
        stop = time.perf_counter()
        complete_elapsed += stop - start 
                    
        accuracy = float(tp + tn)/float(tp + tn + fp + fn)
        precision = float(tp)/float(tp + fp)
        recall = float(tp)/float(tp + fn)
        specificity = float(tn)/float(tn + fp)
        f1score = float(2*tp)/float(2*tp + fp + fn)
        mcc = float(tp*tn - fp*fn)/math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
                
        resultfile.write("Accuracy\t: " + str(accuracy) + "\n")
        resultfile.write("Precision\t: " + str(precision) + "\n")
        resultfile.write("Recall\t\t: " + str(recall) + "\n")
        resultfile.write("Specificity\t: " + str(specificity) + "\n")
        resultfile.write("F1-Score\t: " + str(f1score) + "\n")
        resultfile.write("MCC\t: " + str(mcc) + "\n")
        resultfile.write("Throughput (reads/sec)\t: " + str(float((tp+tn+fp+fn))/complete_elapsed) + "\n")

def compute_minimap_matches_simulated():

    # make one reference file for mapping
    bacreffile = testdir + "ZMCForwardBacReferences.fasta"
    with open(bacreffile, 'wb') as outfile:
        for filename in glob.glob(genomedir + "*.fasta"):
            head, tail = os.path.split(filename)
            if tail.startswith("Saccharomyces"):
                continue
            with open(filename, 'rb') as readfile:
                shutil.copyfileobj(readfile, outfile)

    for sdir in glob.glob(simdir + "ACC*"):
        simdirpath, simdirname = os.path.split(sdir)
        accoutputdir = minioutdir + simdirname + "_360bp/"
        if not os.path.exists(accoutputdir):
	        os.mkdir(accoutputdir)
        tp = 0
        fp = 0
        tn = 0
        fn = 0
        with open(minimapdir +  "minimap2_matching_analysis_" + simdirname + "_360bp.txt", 'w') as resultfile:
            for filename in glob.glob(sdir + "/*.fastq"):
                head, tail = os.path.split(filename)
                newname = os.path.splitext(tail)[0]
                with open(os.devnull, 'w') as devnull:
                    minimapcmd = ["minimap2","-ax","map-ont",bacreffile ,filename]
                    minimapoutput = subprocess.Popen(minimapcmd, stdout=subprocess.PIPE,stderr=devnull)
                    samcmd = ["samtools","view", "-bh","-q","10"]
                    samoutput = subprocess.Popen(samcmd, stdin=minimapoutput.stdout, stdout=subprocess.PIPE, stderr=devnull)
                    samsortcmd = ["samtools", "sort","-O","BAM"]
                    c = subprocess.Popen(samsortcmd ,stdin=samoutput.stdout, stdout=open(accoutputdir + tail + ".bam", 'w'), stderr=devnull)
                    minimapoutput.stdout.close()
                    samoutput.stdout.close()
                    c.communicate()
                    statcommand = ["samtools", "flagstat",accoutputdir + tail + ".bam"]
                    statproc = subprocess.Popen(statcommand, stdout=subprocess.PIPE, stderr=devnull, universal_newlines=True)
                    statout, err = statproc.communicate()
                    statlines = statout.split("\n")
                    mapped = 0
                    for line in statlines:
                        sl = str(line).split()
                        if len(sl) < 3:
                            continue
                        if sl[3] == "mapped":
                            mapped = sl[0]
                        
                    grepcmd = ["grep","-c", "^@",filename]
                    grepoutput = subprocess.Popen(grepcmd, stdout=subprocess.PIPE,stderr=devnull, universal_newlines=True)
                    wcout,err = grepoutput.communicate()
                    nrreads = wcout.split()[0]
                
                    resultfile.write(newname + "\tZMCForwardBacReferences\t" + str(mapped) + "\t" + str(nrreads) + "\n")
                    if not newname.startswith("Saccharomyces"):
                        tp += int(mapped)
                        fn += int(nrreads) - int(mapped)
                    else:
                        tn = int(nrreads) - int(mapped)
                        fp = int(mapped)
                    
            accuracy = float(tp + tn)/float(tp + tn + fp + fn)
            precision = float(tp)/float(tp + fp)
            recall = float(tp)/float(tp + fn)
            specificity = float(tn)/float(tn + fp)
            f1score = float(2*tp)/float(2*tp + fp + fn)
            mcc = float(tp*tn - fp*fn)/math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
                
            resultfile.write("Accuracy\t: " + str(accuracy) + "\n")
            resultfile.write("Precision\t: " + str(precision) + "\n")
            resultfile.write("Recall\t\t: " + str(recall) + "\n")
            resultfile.write("Specificity\t: " + str(specificity) + "\n")
            resultfile.write("F1-Score\t: " + str(f1score) + "\n")
            resultfile.write("MCC\t: " + str(mcc) + "\n")        
    
def main():
    #create_ZMC360bp()
    prepare_seqs_for_indexing()
    create_spumoni_indexes()
    compute_spumoni_pmls_realZMC()
    compute_minimap_matches_realZMC()
    compute_spumoni_pmls_simulated()
    compute_minimap_matches_simulated()

if __name__ == "__main__":
    main()

    
    
    
    
