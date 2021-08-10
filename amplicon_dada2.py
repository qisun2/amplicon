#!/usr/bin/env python3

import logging
from os import listdir
from os.path import isfile, join
import argparse
from pathlib import Path
import os
import sys
from collections import Counter 
from math import log
import subprocess
from operator import itemgetter
import multiprocessing
import re
import traceback
import glob
import pandas as pd
import csv

from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import SeqIO

from collections import defaultdict

#import uuid
#from scipy import stats
#import numpy as np
#import itertools



## dependency on BioHPC
#export PATH=/programs/muscle:$PATH:/programs/bbmap-38.45:$PATH



cutadaptCMD = "cutadapt"
dadaRscript = "dada2.R"


def main():

    commandLine = " ".join(sys.argv[1:])

    global args
    global sampleList
    global sampleToFileList
    global markerList
    global primerFile
    global tagBySampleDir
    global readCountMatrixFile


    parser = argparse.ArgumentParser(description='Run GATK Haplotype Caller on ampseq data.')

    # Required arguments
    parser.add_argument('-s','--sample',type=str,required=True,help='Sample file. Tab delimited text file with 3 or 4 columns: sample_Name, plate_well, fastq_file1, (optional)fastq_file2. Plate_well is a string to uniquely define sample if sample names are duplicated.')
    parser.add_argument('-k','--key',type=str,required=True,help='Key file. Tab delimited text file with 3 columns: amplicon name, 5 prime sequence, 3 prime sequence')
    parser.add_argument('-o','--output',type=str,required=True,help='Output directory')


    # Optional arguments
    parser.add_argument('-i','--skip',type=str,required=False,default="",help='Skip steps. e.g. "-i 1" to skip steps 1. the steps are: 1. split reads by primers; 2. run dada2.R script; 3. final process')
    parser.add_argument('-j','--job',type=int,required=False,default=8,help='Number of simultaneous jobs. Default:8')
    parser.add_argument('-t','--thread',type=int,required=False,default=1,help='Number of threads per job. Default:1')
    parser.add_argument('-l','--minHaplotypeLength',type=int,required=False,default=50,help='Minimum haplotype length (after removing the primers. It must be an integer 1 or larger.) Default:50')
    parser.add_argument('-d','--mergeDuplicate',type=int,required=False,default=0,help='Whether to merge samples with same name. If two samples with same name and same plate_well, they will always be merged. 1: merge; 0: not merge and the duplicated sample will be named <sampleName>__<index starting from 1> . Default:0')
    parser.add_argument('-z','--primerErrorRate',type=float,required=False,default=0.1,help='Mismatch rate between pcr primer and reads, default 0.1')
    parser.add_argument('-r','--refSeq',type=str,required=False,default="",help='Reference sequence for each marker in fasta format')
    parser.add_argument('-e','--evalue',type=float,required=False,default=1e-3,help='Blast to reference sequence evalue cutoff')



    if sys.version_info[0] < 3:
        raise Exception("This code requires Python 3.")

    args=parser.parse_args()

    sampleList = []
    sampleToFileList = []
    markerList= []
    markerDirList = []


    if (not os.path.isfile(args.sample)):
        parser.print_usage()
        print(f"Error: Sample file {args.sample} does not exist!")
        sys.exit()

    if (not os.path.isfile(args.key)):
        parser.print_usage()
        print(f"Error: Key file {args.key} does not exist!")
        sys.exit()





    primerFile = f"{args.output}/primer.fa"
    tagBySampleDir = f"{args.output}/tagBySampleDir"
    readCountMatrixFile = f"{args.output}/markerToSampleReadCountMatrix"



    if (not os.path.exists(args.output)):
        os.mkdir(args.output)
    if (not os.path.exists(tagBySampleDir)):
        os.mkdir(tagBySampleDir)



    logging.basicConfig(filename=args.output + '/run.log',level=logging.DEBUG)

    logging.debug(f"Command: {commandLine}")

    # Process the key file, and prepare linked adapter as required by cutadapt
    primerfh = open(primerFile, "w")
    with open(args.key, 'r') as fhk:
        for line in fhk:
            if (not re.search("\w", line)):
                continue
            line = line.rstrip()
            fieldArray = line.split(sep="\t")
            if (len(fieldArray) < 3):
                print(f"Error: some row in the key file {args.key} have less than 3 columns!")
                sys.exit()
            markerName = re.sub("\W", "", fieldArray[0])
            markerList.append(markerName);
            primer1 = re.sub("\W", "", fieldArray[1])
            primer2 = re.sub("\W", "", fieldArray[2])
            primer2 = revcom(primer2);
            primerfh.write(f">{markerName}\n^{primer1}...{primer2}\n")

            markerDir = f"{args.output}/{markerName}"
            if (not os.path.exists(markerDir)):
                os.mkdir(markerDir)
        fhk.close()
    primerfh.close()

    # check reference sequence file exist if parameter set
    if (args.refSeq != ""):
        if (not os.path.isfile(args.refSeq)):
            parser.print_usage()
            print(f"Error: Reference sequence fasta file {args.refSeq} does not exist!")
            sys.exit()  

    #process sample file, currently, only paired end reads are supported
    ## first check and merge duplicate samples
    checkSampleDup = {}
    fileMerged = {}
    with open(args.sample, 'r') as fhs:
        for line in fhs:
            if (not re.search("\w", line)):
                continue
            line = line.rstrip()
            fieldArray = line.split(sep="\t")
            sampleName = re.sub("\s", "", fieldArray[0])
            plateWell = re.sub("\s", "", fieldArray[1])
            if (args.mergeDuplicate == 0):
                sampleName = sampleName + "__" + plateWell
            if (sampleName in checkSampleDup):
                checkSampleDup[sampleName].append((fieldArray[2], fieldArray[3]))
            else:
                checkSampleDup[sampleName] = [(fieldArray[2], fieldArray[3])]
        fhs.close()

    # execute file merging only if step 1 not skipped

    for sampleName, files in checkSampleDup.items():
        if (len(files)>1):
            ### merge
            file1List = ""
            file2List = ""
            for fff in files:
                file1List+= " " + fff[0]
                file2List+= " " + fff[1]
            
            attachGZ = ""
            if (".gz" in file1List):
                attachGZ = ".gz"

            ## execute merging only if 1 not skipped
            if ("1" not in args.skip):
                mergeCmd1 = f"cat {file1List} > {sampleName}.merged.R1.fastq{attachGZ}"
                mergeCmd2 = f"cat {file2List} > {sampleName}.merged.R2.fastq{attachGZ}"

                logging.debug(f"Merging: {mergeCmd1}")
                logging.debug(f"Merging: {mergeCmd2}")
                os.system(mergeCmd1)
                os.system(mergeCmd2)

            fileMerged[sampleName] = [f"{sampleName}.merged.R1.fastq{attachGZ}", f"{sampleName}.merged.R2.fastq{attachGZ}"]


    with open(args.sample, 'r') as fhs:
        for line in fhs:
            if (not re.search("\w", line)):
                continue
            line = line.rstrip()
            fieldArray = line.split(sep="\t")
            if (len(fieldArray) < 4):
                print(f"Error: Sample file must have at least four columns, and with no header line. Single-end reads are not supported now. Will be added later")
                sys.exit()

            sampleName = re.sub("\s", "", fieldArray[0])
            plateWell = re.sub("\s", "", fieldArray[1])

            if (args.mergeDuplicate == 0):
                sampleName = sampleName + "__" + plateWell


            if (sampleName in fileMerged):
                if (fileMerged[sampleName] == "done"):
                    continue
                else:
                    fieldArray[2] = fileMerged[sampleName][0]
                    fieldArray[3] = fileMerged[sampleName][1]
                    fileMerged[sampleName] = "done"

            if (sampleName in sampleList):
                continue
            sampleList.append(sampleName)

            if ((not os.path.isfile(fieldArray[2])) and  ("1" not in args.skip)):
                print(f"Error: Sample fastq file {fieldArray[2]} does not exist!")
                sys.exit()
            if ((not os.path.isfile(fieldArray[3])) and ("1" not in args.skip)):
                print(f"Error: Sample fastq file {fieldArray[3]} does not exist!")
                sys.exit()
                
            sampleToFileList.append((sampleName, fieldArray[2], fieldArray[3]))

    # split the read by primers, and remove primer from each read, and collapase identical reads, and keep top <arg.maxHaplotypePerSample> tags per sample
    if ("1" not in args.skip):
        logging.info("Step 1: splitByPrimer")
        pool = multiprocessing.Pool(processes= args.job)
        pool.starmap(splitByCutadapt, sampleToFileList, chunksize=20)
        pool.close()

    # get sequence tag list across population
    if ("2" not in args.skip):
        logging.info("Step 2: run dada2")
        for markerName in markerList:
            runDada(markerName)


    if ("3" not in args.skip):
        logging.info("Step 3: final process")

        # process reference sequence file
        if (args.refSeq != ""):            
            refseqExist = set()
            for record in SeqIO.parse(args.refSeq, "fasta"):
                markerName = record.id

                if (markerName in markerList):
                    refseqExist.add(markerName)
                    wsh = open(f"{args.output}/{markerName}/refseq.fas", "w")
                    SeqIO.write(record, wsh, "fasta")
                    wsh.close()
                    cmd = f"makeblastdb -in {args.output}/{markerName}/refseq.fas -dbtype nucl"
                    returned_value = subprocess.call(cmd, shell=True)

            for m in markerList:
                if m not in refseqExist:
                    print(f"Error: Reference sequence fasta file {args.refSeq} does not contain reference sequence for marker {m}!")
                    sys.exit()
        
        markerDirList = []
        for markerName in markerList:
            markerDirList.append((markerName,))

        pool = multiprocessing.Pool(processes= args.job)
        
        pool.starmap(finalProcess, markerDirList)
        pool.close()


def splitByCutadapt(sampleName, file1, file2):
    try:
        sampleDir = f"{args.output}/{sampleName}"
        if (not os.path.exists(sampleDir)):
            os.mkdir(sampleDir)

        #run cutadapt to demultiplexing by primers
        cmd = f"{cutadaptCMD} --quiet -e {args.primerErrorRate} --minimum-length={args.minHaplotypeLength} --discard-untrimmed -g file:{primerFile} -o {sampleDir}/{{name}}_R1.fastq -p {sampleDir}/{{name}}_R2.fastq {file1} {file2} "
        returned_value = subprocess.call(cmd, shell=True)
        logging.info(f"demultiplexing {sampleName} done: {returned_value}")


        # move splitted file to corresponding marker directory
        for markerName in markerList:
            if (os.path.exists(f"{sampleDir}/{markerName}_R1.fastq")):
                cmd1 = f"mv {sampleDir}/{markerName}_R1.fastq {args.output}/{markerName}/{sampleName}_R1.fastq"
                cmd2 = f"mv {sampleDir}/{markerName}_R2.fastq {args.output}/{markerName}/{sampleName}_R2.fastq"

                subprocess.call(cmd1, shell=True)
                subprocess.call(cmd2, shell=True)


        os.system("rm -rf %s" % sampleDir)
        return 1 
    except Exception as e:
        print(f'Caught exception in job {sampleName} ')
        traceback.print_exc()
        print()
        raise e


def runDada(markerPath):
    try:
        cmd = f"{dadaRscript} {args.output}/{markerPath}"
        returned_value = subprocess.call(cmd, shell=True)
        return 1 
    except Exception as e:
        print(f'Caught exception in job {markerPath} ')
        traceback.print_exc()
        print()
        raise e


def finalProcess(markerPath):
    try:
        if (args.refSeq != ""):
            acceptedHaplotype= set()
            acceptedHaplotypeSeqs= set()

            fastaFile = f"{args.output}/{markerPath}/{markerPath}.uniqueSeqs.fasta"
            fastaFilteredFile = f"{args.output}/{markerPath}/{markerPath}.uniqueSeqs.filtered.fasta"
            blastResult = f"{args.output}/{markerPath}/blastresults"


            if (os.path.isfile(fastaFile)):
                cmd = f"blastn -query {fastaFile} -db {args.output}/{markerPath}/refseq.fas  -dust no -outfmt \"6 std qcovs\"  -out {blastResult}"
                returned_value = subprocess.call(cmd, shell=True)

                with open (blastResult, "r") as BFH:
                    for line in BFH:
                        line = line.rstrip()
                        qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qcov= line.split(sep="\t")
                        evalue =float(evalue)
                        if (evalue< args.evalue):
                            acceptedHaplotype.add(qseqid)
                        
                acceptedRecords = []
                for record in SeqIO.parse(fastaFile, "fasta"):
                    if (record.id in acceptedHaplotype):
                        acceptedHaplotypeSeqs.add(str(record.seq))
                        acceptedRecords.append(record)
                wsh = open(fastaFilteredFile, "w")
                SeqIO.write(acceptedRecords, wsh, "fasta")
                wsh.close()

        ## filter .seqtab.csv file
        csvFile = f"{args.output}/{markerPath}/{markerPath}.seqtab.csv"
        modCsvFile = f"{args.output}/{markerPath}/{markerPath}.seqtab.mod.csv"

        dataMatrix= pd.read_csv(csvFile, header=0, index_col=0)

        #modify column index
        samplesInFile = [item.replace("_F_filt.fastq.gz","") for item in dataMatrix.columns]
        dataMatrix.columns = samplesInFile

        dataMatrix = dataMatrix.loc[acceptedHaplotypeSeqs]

        dataMatrix = dataMatrix.reindex(columns=sampleList, fill_value=0)
        
        dataMatrix.to_csv(modCsvFile, quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
        
        return 1 
    except Exception as e:
        print(f'Caught exception in job {markerPath} ')
        traceback.print_exc()
        print()
        raise e


def revcom(inputSeq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'} 
    seq = re.sub("\s", '', inputSeq.upper())
    seq = re.sub("[^ACGT]", 'N', seq)
    bases = list(seq[::-1]) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)

def checkApp(checkCommand):
    p = subprocess.Popen(f"which {checkCommand}", stdout=subprocess.PIPE, shell=True)
    res = p.stdout.readlines()
    p.stdout.close()

    if (len(res)> 0): 
        print (f"{checkCommand}: YES")
        return True
    else:
        print (f"{checkCommand}: NOT Available")
        return False




if __name__=="__main__":
    main()


