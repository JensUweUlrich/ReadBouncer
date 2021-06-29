#!/usr/bin/python3
import shutil
import glob
import os
import subprocess
import math

rbpath = "/usr/local/ReadBouncer/bin/ReadBouncer"
testdir = "testData/"
genomedir = "Genomes/"
simdir = testdir + "SimulatedReads/"
zmc360 = testdir + "RealMock_360bp/"

resultdir = "results/"
rbdir = resultdir + "ReadBouncer/"


if not os.path.exists(rbdir):
	os.mkdir(rbdir)

def buildIBF():
    buildcmd = [rbpath, "ibfbuild", "-i", testdir + "ZMCForwardBacReferences.fasta", "-o", testdir + "ZMCForwardBacReferences.ibf"]
    subprocess.run(buildcmd)


def classifyRealZMCreads():
    tp = 0
    fp = 0
    tn = 0
    fn = 0
    all_reads = 0
    peakrss = 0
    with open(rbdir +  "ReadBouncer_classify_analysis_ZMCReal_360bp.txt", 'w') as resultfile:
        complete_elapsed = 0.0
        
        for	filename in glob.glob(zmc360 + "*.fasta"):
            head, tail = os.path.split(filename)
            newname = os.path.splitext(tail)[0]

            classifycmd = [rbpath, "classify", "-d", testdir + "ZMCForwardBacReferences.ibf", "-r", filename, "-p", "360"]
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
                if sl[0] == "Peak":
                    if int(sl[3]) > peakrss:
                        peakrss = int(sl[3])

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
        resultfile.write("Throughput in reads/sec\t: " + str(float(all_reads)/complete_elapsed) + "\n")
        resultfile.write("Peak RSS in MegaByte\t: " + str(peakrss) + "\n")


def classifySimulatedReads():

    for sdir in glob.glob(simdir + "ACC*"):
        simdirpath, simdirname = os.path.split(sdir)
        tp = 0
        fp = 0
        tn = 0
        fn = 0
        with open(rbdir +  "ReadBouncer_classify_analysis_" + simdirname + "_360bp.txt", 'w') as resultfile:
            for filename in glob.glob(sdir + "/*.fastq"):
                head, tail = os.path.split(filename)
                newname = os.path.splitext(tail)[0]

                classifycmd = [rbpath, "classify", "-d", testdir +  "ZMCForwardBacReferences.ibf", "-r", filename]
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
    buildIBF()
    classifyRealZMCreads()
    classifySimulatedReads()

if __name__ == "__main__":
    main()
