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
# download reference sequence from NCBI and put in genome dir
mouse_ref = genomedir + "GCF_000001635.27_GRCm39_genomic.fna"
spumoni_posIndex = "mouse_positive_index.fasta"
spumoni_nullIndex = "mouse_null_index.fasta"

rattus_reads = "testData/rattus_360bp_reads.fastq"
# copy from ZymoMockCommunity test data
zymo_reads = "testData/ZymoMCreads.fastq"

def simulate_rat_reads():
    pbsimcmd = ["pbsim", "--depth", "2", "--accuracy-mean", "0.9","--length-min","360", "--length-max", "360", "--hmm_model", pbmodeldir + "R94.model", "--seed" , "20210617", filename]
    subprocess.run(pbsimcmd)
    with open(rattus_reads, 'wb') as outfile:
        for filename in glob.glob("sd_*.fastq"):
            if filename == outfilename:
                # don't want to copy the output into the output
                continue
            with open(filename, 'rb') as readfile:
                shutil.copyfileobj(readfile, outfile)
        
    for file in glob.glob("sd_*"):
       	os.remove(file)

# create human references with their reverse complement sequences
# add forward and reverse complement of every reference sequence to the files
# spumoni needs both for index building
def prepare_seqs_for_indexing():

    # make one reference file for mapping
    outfilename = testdir + spumoni_posIndex
    with open(outfilename, 'w') as f_out:
        for record in SeqIO.parse(mouse_ref, "fasta"):
            new_record = SeqRecord(seq=record.seq.back_transcribe(), id=record.id)
            r=SeqIO.write(new_record, f_out, 'fasta')
            if r!=1: print('Error while writing sequence:  ' + record.id)    
            new_record = SeqRecord(seq=record.seq.back_transcribe().reverse_complement(), id=record.id + "_revComp")
            r=SeqIO.write(new_record, f_out, 'fasta')
            if r!=1: print('Error while writing sequence:  ' + new_record.id)
     

def create_spumoni_indexes():
    #build positive index for references
    monibuildcmd = ["python3",monipath, "build", "-r", testdir + spumoni_posIndex, "-f", "--spumoni", "-t", "8"]
    subprocess.run(monibuildcmd)
    for file in glob.glob(testdir + spumoni_posIndex + ".*"):
        shutil.move(file, indexdir)

    # build null index for references
    with open(testdir + spumoni_nullIndex, 'w') as f_out:
        for record in SeqIO.parse(testdir + spumoni_posIndex, "fasta"):
            new_record = SeqRecord(seq=record.seq[::-1], id=record.id)
            r=SeqIO.write(new_record, f_out, 'fasta')
            if r!=1: print('Error while writing sequence:  ' + record.id)

    monibuildcmd = ["python3",monipath, "build", "-r", testdir + spumoni_nullIndex, "-f", "--spumoni","-t", "8"]
    subprocess.run(monibuildcmd)
    for file in glob.glob(testdir + spumoni_nullIndex + ".*"):
        shutil.move(file, indexdir)
        

def compute_spumoni_pmls():

    tp = 0
    fp = 0
    tn = 0
    fn = 0
    complete_elapsed = 0.0
    org_counter = 0
    nr_allreads = 0
    with open(spumonidir +  "Spumoni_pml_analysis.txt", 'w') as resultfile:
        head, tail = os.path.split(rattus_reads)
        newname = os.path.splitext(tail)[0]
        org_counter += 1
        start = time.perf_counter()
        # pseudo-ms for positive index of bacterial references
        pseudomscmd = ["python3", monipath, "pseudo-ms", "-i", indexdir + spumoni_posIndex, "-p", rattus_reads]
        subprocess.run(pseudomscmd)
        # pseudo-ms for null index of bacterial references
        pseudomscmd = ["python3", monipath, "pseudo-ms", "-i", indexdir + spumoni_nullIndex, "-p", rattus_reads]
        subprocess.run(pseudomscmd)
        for file in glob.glob(rattus_reads + "_*"):
            shutil.move(file, outdir)
        os.remove(rattus_reads + ".pseudo_ms.log")           
                        
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
                        
        resultfile.write(newname + "\trattus_norvegicus\t" + str(nrmatches) + "\t" + str(nrreads) + "\n")
        nr_allreads += nrreads
        tp += int(nrmatches)
        fn += int(nrreads) - int(nrmatches)

        head, tail = os.path.split(zymo_reads)
        newname = os.path.splitext(tail)[0]    
        # pseudo-ms for positive index of yeast
        pseudomscmd = ["python3", monipath, "pseudo-ms", "-i", indexdir + spumoni_posIndex, "-p", zymo_reads]
        subprocess.run(pseudomscmd)
        # pseudo-ms for null index of yeast references
        pseudomscmd = ["python3", monipath, "pseudo-ms", "-i", indexdir + spumoni_nullIndex, "-p", zymo_reads]
        subprocess.run(pseudomscmd)
        for file in glob.glob(zymo_reads + "_*"):
            shutil.move(file, outdir)
        os.remove(zymo_reads + ".pseudo_ms.log")  
                
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
                        
        resultfile.write(newname + "\tZymoMC\t" + str(nrmatches) + "\t" + str(nrreads) + "\n")
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
    mapper = mappy.Aligner(mouse_ref, preset="map-ont", n_threads=1)

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
        head, tail = os.path.split(rattus_reads)
        newname = os.path.splitext(tail)[0]
        for record in SeqIO.parse(rattus_reads, "fastq"):
            nrreads += 1
            for hit in mapper.map(record.seq):
                mapped += 1
                break
        
        resultfile.write(newname + "\trattus_norvegicus\t" + str(mapped) + "\t" + str(nrreads) + "\n")
        tp += int(mapped)
        fn += int(nrreads) - int(mapped)
        nr_allreads = nrreads
        
        # map microbiome reads 
        
        mapped = 0
        nrreads = 0
        head, tail = os.path.split(zymo_reads)
        newname = os.path.splitext(tail)[0]
        for record in SeqIO.parse(zymo_reads, "fastq"):
            nrreads += 1
            for hit in mapper.map(record.seq):
                mapped += 1
                break
        
        stop = time.perf_counter()
        complete_elapsed += stop - start    
        resultfile.write(newname + "\tZymoMC\t" + str(mapped) + "\t" + str(nrreads) + "\n")
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
    buildcmd = [rbpath, "ibfbuild", "-i", mouse_ref, "-o", testdir + "MusMusculus.ibf"]
    subprocess.run(buildcmd)
    
def classifyReadBouncerReads():
    tp = 0
    fp = 0
    tn = 0
    fn = 0
    all_reads = 0
    with open(rbdir +  "ReadBouncer_classify_analysis.txt", 'w') as resultfile:
        complete_elapsed = 0.0
        
        mapped = 0
        nrreads = 0
        head, tail = os.path.split(rattus_reads)
        newname = os.path.splitext(tail)[0]
        
        complete_elapsed_2 = 0.0
        start = time.perf_counter()
        classifycmd = [rbpath, "classify", "-d", testdir + "MusMusculus.ibf", "-r", rattus_reads]
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

        resultfile.write(newname + "\trattus_norvegicus\t" + str(mapped) + "\t" + str(nrreads) + "\n")
        tp += int(mapped)
        fn += int(nrreads) - int(mapped)
        
        mapped = 0
        nrreads = 0
        head, tail = os.path.split(zymo_reads)
        newname = os.path.splitext(tail)[0]

        classifycmd = [rbpath, "classify", "-d", testdir + "MusMusculus.ibf", "-r", zymo_reads]
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

        resultfile.write(newname + "\tZymoMC\t" + str(mapped) + "\t" + str(nrreads) + "\n")
        fp += int(mapped)
        tn += int(nrreads) - int(mapped)
        stop = time.perf_counter()
        complete_elapsed_2 += stop - start  
            
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
        resultfile.write("Throughput in reads/sec\t: " + str(float(all_reads)/complete_elapsed) + "\n")
        resultfile.write("Throughput in reads/sec\t: " + str(float(all_reads)/complete_elapsed_2) + "\n")


    
def main():
    simulate_rat_reads()
    prepare_seqs_for_indexing()
    create_spumoni_indexes()
    compute_spumoni_pmls()
    compute_minimap_matches()
    buildIBF()
    classifyNanoLiveReads()
    

if __name__ == "__main__":
    main()

    
    
    
    
