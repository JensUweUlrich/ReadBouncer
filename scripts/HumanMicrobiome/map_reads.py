import shutil
import glob
import os
import subprocess
import time
import math
import mappy
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

testdir = "testData/"
genomedir = "Genomes/"

# paths have to be set accordingly before running
monipath = "/home/jens/tools/spumoni/build/moni"
analyzepath = "/home/jens/tools/spumoni/analysis/analyze_pml.py"
rbpath = "/usr/local/ReadBouncer/bin/ReadBouncer"

resultdir = "results/"
spumonidir = resultdir + "spumoni/"
indexdir = "index/"
outdir = spumonidir + "output/"

minimapdir = resultdir + "minimap2/"
minioutdir = minimapdir + "output/"
rbdir = resultdir + "ReadBouncer/"


if not os.path.exists(resultdir):
	os.mkdir(resultdir)
if not os.path.exists(spumonidir):
	os.mkdir(spumonidir)
if not os.path.exists(indexdir):
	os.mkdir(indexdir)
if not os.path.exists(outdir):
	os.mkdir(outdir)
if not os.path.exists(minimapdir):
	os.mkdir(minimapdir)
if not os.path.exists(minioutdir):
	os.mkdir(minioutdir)
if not os.path.exists(rbdir):
	os.mkdir(rbdir)

# reference sequences for depletion
chm13ref = "chm13.draft_v1.0.fasta"
# download from https://genome-idx.s3.amazonaws.com/spumoni/one_human_genome/one_human_genome_positive_index.spumoni
spumoni_posIndex = "chm13_positive_index.fasta"
# download from https://genome-idx.s3.amazonaws.com/spumoni/one_human_genome/one_human_genome_null_index.spumoni
spumoni_nullIndex = "chm13_null_index.fasta"

# have to be created in testData dir by using create360bpReads.py on both read data sets
human_reads = testdir + "HumanGenomeReads_360bp.fasta"
mircobiome_reads = testdir + "CleanMicrobiome_360bp.fasta" 

def compute_spumoni_pmls():

    tp = 0
    fp = 0
    tn = 0
    fn = 0
    complete_elapsed = 0.0
    org_counter = 0
    nr_allreads = 0
    with open(spumonidir +  "pml_analysis_HumanReads_360bp.txt", 'w') as resultfile:
        head, tail = os.path.split(human_reads)
        newname = os.path.splitext(tail)[0]
        org_counter += 1
        start = time.perf_counter()
        # pseudo-ms for positive index of bacterial references
        pseudomscmd = ["python3", monipath, "pseudo-ms", "-i", indexdir + spumoni_posIndex, "-p", human_reads]
        subprocess.run(pseudomscmd)
        # pseudo-ms for null index of bacterial references
        pseudomscmd = ["python3", monipath, "pseudo-ms", "-i", indexdir + spumoni_nullIndex, "-p", human_reads]
        subprocess.run(pseudomscmd)
        for file in glob.glob(human_reads + "_*"):
            shutil.move(file, outdir)
        os.remove(human_reads + ".pseudo_ms.log")           
                        
        analyzecmd = ["python3", analyzepath, "-p", outdir + tail + "_" + spumoni_posIndex + ".pseudo_lengths", "-n", outdir + tail + "_" + spumoni_nullIndex + ".pseudo_lengths", "-r", "90"]
        proc = subprocess.Popen(analyzecmd , stdout=open(outdir + newname + "_pml_analysis.txt", 'w'))
        proc.communicate()    

        stop = time.perf_counter()
        complete_elapsed += stop - start    
        nrmatches = 0
        nrreads = 0
        with open(outdir + newname + "_pml_analysis.txt", 'r') as f:
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
                        
        resultfile.write(newname + "\tchm13_v1\t" + str(nrmatches) + "\t" + str(nrreads) + "\n")
        nr_allreads += nrreads
        tp += int(nrmatches)
        fn += int(nrreads) - int(nrmatches)

        head, tail = os.path.split(mircobiome_reads)
        newname = os.path.splitext(tail)[0]    
        # pseudo-ms for positive index chm13
        pseudomscmd = ["python3", monipath, "pseudo-ms", "-i", indexdir + spumoni_posIndex, "-p", mircobiome_reads]
        subprocess.run(pseudomscmd)
        # pseudo-ms for null index of chm13
        pseudomscmd = ["python3", monipath, "pseudo-ms", "-i", indexdir + spumoni_nullIndex, "-p", mircobiome_reads]
        subprocess.run(pseudomscmd)
        for file in glob.glob(mircobiome_reads + "_*"):
            shutil.move(file, outdir)
        os.remove(mircobiome_reads + ".pseudo_ms.log")  
                
        analyzecmd = ["python3", analyzepath, "-p", outdir + tail + "_" + spumoni_posIndex + ".pseudo_lengths", "-n", outdir + tail + "_" + spumoni_nullIndex + ".pseudo_lengths", "-r", "90"]
        proc = subprocess.Popen(analyzecmd , stdout=open(outdir + newname + "_pml_analysis.txt", 'w'))
        proc.communicate()
        nrmatches = 0
        nrreads = 0
        with open(outdir + newname + "_pml_analysis.txt", 'r') as f:
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
                        
        resultfile.write(newname + "\tchm13_v1\t" + str(nrmatches) + "\t" + str(nrreads) + "\n")
        nr_allreads += nrreads
        fp += int(nrmatches)
        tn += int(nrreads) - int(nrmatches)
            
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
        resultfile.write("Throughput (reads/sec)\t: " + str(float(nr_allreads)/complete_elapsed) + "\n")
    

def compute_minimap_matches():

    # create mapper object
    mapper = mappy.Aligner(chm13ref, preset="map-ont",n_threads=1)

    tp = 0
    fp = 0
    tn = 0
    fn = 0
    complete_elapsed = 0.0
    with open(minimapdir +  "minimap2_matching_analysis.txt", 'w') as resultfile:

        start = time.perf_counter()
        # map human reads
        mapped = 0
        nrreads = 0
        head, tail = os.path.split(human_reads)
        newname = os.path.splitext(tail)[0]
        for record in SeqIO.parse(human_reads, "fasta"):
            nrreads += 1
            for hit in mapper.map(record.seq):
                mapped += 1
                break
        
        resultfile.write(newname + "\tchm13_v1\t" + str(mapped) + "\t" + str(nrreads) + "\n")
        tp += int(mapped)
        fn += int(nrreads) - int(mapped)
        nr_allreads = nrreads
        
        # map microbiome reads 
        
        mapped = 0
        nrreads = 0
        head, tail = os.path.split(mircobiome_reads)
        newname = os.path.splitext(tail)[0]
        for record in SeqIO.parse(mircobiome_reads, "fasta"):
            nrreads += 1
            for hit in mapper.map(record.seq):
                mapped += 1
                break
        
        stop = time.perf_counter()
        complete_elapsed += stop - start    
        resultfile.write(newname + "\tchm13_v1\t" + str(mapped) + "\t" + str(nrreads) + "\n")
        fp += int(mapped)
        tn += int(nrreads) - int(mapped)
        nr_allreads += nrreads
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
        resultfile.write("Throughput (reads/sec)\t: " + str(float(nr_allreads)/complete_elapsed) + "\n")

def buildIBF():
    buildcmd = [rbpath, "ibfbuild", "-i", chm13ref, "-o", testdir + "chm13.ibf"]
    subprocess.run(buildcmd)
    
def classifyReadBouncerReads():
    tp = 0
    fp = 0
    tn = 0
    fn = 0
    all_reads = 0
    with open(rbdir +  "ReadBouncer_classify_analysis_360bp.txt", 'w') as resultfile:
        complete_elapsed = 0.0
        
        mapped = 0
        nrreads = 0
        head, tail = os.path.split(human_reads)
        newname = os.path.splitext(tail)[0]

        classifycmd = [rbpath, "classify", "-d", testdir + "chm13.ibf", "-r", human_reads]
        classifyproc = subprocess.Popen(classifycmd, stdout=subprocess.PIPE, universal_newlines=True)
        out, err = classifyproc.communicate()
        statlines = out.split("\n")
        mapped = 0
        nrreads = 0
        for line in statlines:
            sl = str(line).split()
            if len(sl) < 2:
                continue
            if sl[2] == "classified":
                mapped = sl[5]
            if sl[2] == "all":
                nrreads = sl[5]
                all_reads += int(nrreads)
            if sl[0] == "Real":
                complete_elapsed += float(sl[3])

        resultfile.write(newname + "\tchm13_v1\t" + str(mapped) + "\t" + str(nrreads) + "\n")
        tp += int(mapped)
        fn += int(nrreads) - int(mapped)
        
        mapped = 0
        nrreads = 0
        head, tail = os.path.split(mircobiome_reads)
        newname = os.path.splitext(tail)[0]

        classifycmd = [rbpath, "classify", "-d", testdir + "chm13.ibf", "-r", mircobiome_reads]
        classifyproc = subprocess.Popen(classifycmd, stdout=subprocess.PIPE, universal_newlines=True)
        out, err = classifyproc.communicate()
        statlines = out.split("\n")
        mapped = 0
        nrreads = 0
        for line in statlines:
            sl = str(line).split()
            if len(sl) < 2:
                continue
            if sl[2] == "classified":
                mapped = sl[5]
            if sl[2] == "all":
                nrreads = sl[5]
                all_reads += int(nrreads)
            if sl[0] == "Real":
                complete_elapsed += float(sl[3])

        resultfile.write(newname + "\tchm13_v1\t" + str(mapped) + "\t" + str(nrreads) + "\n")
        fp += int(mapped)
        tn += int(nrreads) - int(mapped)
        
            
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
        resultfile.write("Throughput in bp/sec\t: " + str(float(all_reads)/complete_elapsed) + "\n")


    
def main():
    compute_spumoni_pmls()
    compute_minimap_matches()
    buildIBF()
    classifyNanoLiveReads()
    

if __name__ == "__main__":
    main()

    
    
    
    
