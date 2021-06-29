#!/usr/bin/python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import getopt,sys,re
from csv import reader
#import pysam
import numpy as np
#import matplotlib.pyplot as plt
#import pandas as pd
#import seaborn as sns

def usage():
    print("Usage: getRunStats.py -s <ReadFish decision stats> -c <ReadBouncer decision stats>")

def quantile_normalize(df):
    """
    input: dataframe with numerical columns
    output: dataframe with quantile normalized values
    """
    df_sorted = pd.DataFrame(np.sort(df.values,
                                     axis=0), 
                             index=df.index, 
                             columns=df.columns)
    df_mean = df_sorted.mean(axis=1)
    df_mean.index = np.arange(1, len(df_mean) + 1)
    df_qn =df.rank(method="min").stack().astype(int).map(df_mean).unstack()
    return(df_qn)

try:
    options, remainder=getopt.getopt(sys.argv[1:], 'c:s:h')

except getopt.GetoptError as err:
    print(str(err))
    usage()
    sys.exit()

for opt, arg in options:
    if opt in ('-s'):
        readfish_file=arg
    if opt in ('-h'):
        usage()
        sys.exit()
    elif opt in ('-c'):
        readbouncer_file=arg

readfish = {}
duplications = 0
unblock_reads = 0
avg_unblock_len = 0.0
stop_reads = 0
avg_stop_len = 0.0
readfishlen_arr = []
with open(readfish_file, 'r') as csv_in:
    csv_reader = reader(csv_in, delimiter='\t')
    header = next(csv_reader)
    print(header)
   
    for row in csv_reader:
        read_id = row[2]#row[2] #row[5] #row[1]
        decision = row[8]#row[8] #row[7] #row[5]
        seqlen = row[5]#row[5] #row[6] #row[4]
        if read_id in readfish.keys():
            duplications += 1
            continue
        if decision == "unblock":
            unblock_reads += 1
            readfish[read_id] = {}
            readfish[read_id]["unblock"] = True
            readfish[read_id]["RUlen"] = seqlen
            readfishlen_arr.append(float(seqlen))
            #avg_unblock_len += float(seqlen)
        elif decision == "exceeded_max_chunks_unblocked":
            stop_reads += 1
            readfish[read_id] = {}
            readfish[read_id]["unblock"] = False
#            reads[read_id]["RUlen"] = seqlen
#            avg_stop_len += float(seqlen)
        #elif decision != "no_decision":
            #print(decision)
#avg_unblock_len /= float(unblock_reads)
np_arr = np.array(readfishlen_arr)
#avg_stop_len /= float(stop_reads)
print("Number of duplications (Readfish)            : " + str(duplications))
print("Number of unblocked reads (Readfish)         : " + str(unblock_reads))
print("Average read length (Readfish)               : " + str(np.mean(np_arr)))
print("Std. Dev. read length (Readfish)             : " + str(np.std(np_arr)))
print("Median read length (Readfish)                : " + str(np.median(np_arr)))
print("Max read length (Readfish)                   : " + str(np.max(np_arr)))
print("Min read length (Readfish)                   : " + str(np.min(np_arr)))
print("Number unclassified (Readfish)               : " + str(stop_reads))
readbouncer = {}
duplications = 0
unblock_reads = 0
avg_unblock_len = 0.0
stop_reads = 0
avg_stop_len = 0.0
readbouncerlen_arr = []
with open(readbouncer_file, 'r') as csv_in:
    csv_reader = reader(csv_in, delimiter=';')
    header = next(csv_reader)
    print(header)
   
    for row in csv_reader:
        read_id = row[1]#row[2] #row[5] #row[1]
        decision = row[5]#row[8] #row[7] #row[5]
        seqlen = row[4]#row[5] #row[6] #row[4]
        if read_id in readbouncer.keys():
            duplications += 1
            continue
        if decision == "unblock":
            unblock_reads += 1
            readbouncer[read_id] = {}
            readbouncer[read_id]["unblock"] = True
            readbouncer[read_id]["RUlen"] = seqlen
            readbouncerlen_arr.append(float(seqlen))
            #avg_unblock_len += float(seqlen)
        elif decision == "stop_receiving":
            stop_reads += 1
np_arr = np.array(readbouncerlen_arr)
#avg_stop_len /= float(stop_reads)
print("Number of duplications (ReadBouncer)         : " + str(duplications))
print("Number of unblocked reads (ReadBouncer)      : " + str(unblock_reads))
print("Average read length (ReadBouncer)            : " + str(np.mean(np_arr)))
print("Std. Dev. read length (ReadBouncer)          : " + str(np.std(np_arr)))
print("Median read length (ReadBouncer)             : " + str(np.median(np_arr)))
print("Max read length (ReadBouncer)                : " + str(np.max(np_arr)))
print("Min read length (ReadBouncer)                : " + str(np.min(np_arr)))
print("Number unclassified (ReadBouncer)            : " + str(stop_reads))
#plt.hist(np_arr, density=False, bins=320)  # density=False would make counts
#plt.ylabel('Read Counts')
#plt.xlabel('Read length')
#plt.axis(xmin=0,xmax=1600) 
#plt.show()
#df=pd.DataFrame({"Readfish":readfishlen_arr})
#df_qn=quantile_normalize(df)
#both = [readfishlen_arr,
#     readbouncerlen_arr]
#sns.boxplot(data=both)
# set x-axis label
#plt.xlabel("Samples", size=18)
# set y-axis label
#plt.ylabel("Normalized Read Counts", size=18)
#plt.title("Boxplot after Quantile Normalization")
#plt.show()
#plt.savefig('Boxplot_before_Quantile_Normalization_Seaborn.png',dpi=150)

