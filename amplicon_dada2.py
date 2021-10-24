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
import swalign

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
blastCmd = "blastn"
makeblastdbCmd = "makeblastdb"

def main():

    commandLine = " ".join(sys.argv[1:])

    global args
    global sampleList
    global sampleToFileList
    global markerList
    global primerFileF
    global primerFileR
    global tagBySampleDir
    global readCountMatrixFile
    global sw
    global refseq_5p25nt_dict
    global refseq_3p25nt_dict


    parser = argparse.ArgumentParser(description='Run GATK Haplotype Caller on ampseq data.')

    # Required arguments
    parser.add_argument('-s','--sample',type=str,required=True,help='Sample file. Tab delimited text file with 3 or 4 columns: sample_Name, plate_well, fastq_file1, (optional)fastq_file2. Plate_well is a string to uniquely define sample if sample names are duplicated.')
    parser.add_argument('-k','--key',type=str,required=True,help='Key file. Tab delimited text file with 3 columns: amplicon name, 5 prime sequence, 3 prime sequence')
    parser.add_argument('-o','--output',type=str,required=True,help='Output directory')


    # Optional arguments
    parser.add_argument('-i','--skip',type=str,required=False,default="",help='Skip steps. e.g. "-i 1" to skip steps 1. the steps are: 1. split reads by primers; 2. run dada2.R script; 3. final process; 4. delete intermediate files')
    parser.add_argument('-j','--job',type=int,required=False,default=8,help='Number of simultaneous jobs. Default:8')
    parser.add_argument('-t','--thread',type=int,required=False,default=1,help='Number of threads per job. Default:1')
    parser.add_argument('-l','--minHaplotypeLength',type=int,required=False,default=50,help='Minimum haplotype length (after removing the primers. It must be an integer 1 or larger.) Default:50')
    parser.add_argument('-d','--mergeDuplicate',type=int,required=False,default=0,help='Whether to merge samples with same name. If two samples with same name and same plate_well, they will always be merged. 1: merge; 0: not merge and the duplicated sample will be named <sampleName>__<index starting from 1> . Default:0')
    parser.add_argument('-z','--primerErrorRate',type=float,required=False,default=0.15,help='Mismatch rate between pcr primer and reads, default 0.1')
    parser.add_argument('-f','--filter',type=str,required=False,default="1",help='Filter alleles by reference sequence. 1. No filter; 2 Filter by smith-waterman alignment to reference; 3. Filter by BLAST. default 1')
    parser.add_argument('-x','--alnpct',type=float,required=False,default=0.7,help='percent identity of alignment')
    parser.add_argument('-r','--refSeq',type=str,required=False,default="",help='Reference sequence for each marker in fasta format. If the value is set as 1, the first allele from dada2 is used as refseq')
    parser.add_argument('-e','--evalue',type=float,required=False,default=1e-2,help='Blast to reference sequence evalue cutoff')



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





    primerFileF = f"{args.output}/primerF.fa"
    primerFileR = f"{args.output}/primerR.fa"
    tagBySampleDir = f"{args.output}/tagBySampleDir"
    readCountMatrixFile = f"{args.output}/markerToSampleReadCountMatrix"



    if (not os.path.exists(args.output)):
        os.mkdir(args.output)
    if (not os.path.exists(tagBySampleDir)):
        os.mkdir(tagBySampleDir)



    logging.basicConfig(filename=args.output + '/run.log',level=logging.DEBUG)

    logging.debug(f"Command: {commandLine}")

    # Process the key file, and prepare linked adapter as required by cutadapt
    primerfhF = open(primerFileF, "w")
    primerfhR = open(primerFileR, "w")
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
            primerfhF.write(f">{markerName}\n^{primer1}\n")
            primerfhR.write(f">{markerName}\n^{primer2}\n")

            markerDir = f"{args.output}/{markerName}"
            if (not os.path.exists(markerDir)):
                os.mkdir(markerDir)
        fhk.close()
    primerfhF.close()
    primerfhR.close()

    # check reference sequence file exist if parameter set
    if (args.filter not in ["1", "2", "3"]):
        print (f"Error: The --filter option must be a value of 1, 2 or 3")
        parser.print_usage()
        sys.exit()  


    if (args.filter == "1"):
        pass
    else:
        if (args.refSeq == ""):
            print(f"Error: Reference sequence fasta file {args.refSeq} are not set!")
            parser.print_usage()
            sys.exit()  
        elif (args.refSeq != "1") and (not os.path.isfile(args.refSeq)):
            print(f"Error: Reference sequence fasta file {args.refSeq} does not exist!")
            parser.print_usage()
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
            logging.info (f"start dada2 on marker {markerName}")
            returnValue= runDada(markerName)
            logging.info (f"return dada run value for  {markerName} is {returnValue}")
            if (returnValue !=0):
                logging.info (f"Rerun the marker {markerName}")
                runDada(markerName)


    if ("3" not in args.skip):
        logging.info("Step 3: final process")

        # process reference sequence file
        if (args.filter != "1"):
            # create a refseq if refseq input value as 1
            refseqFile = args.refSeq
            if (args.refSeq == "1"):
                refseqFile = f"{args.output}/refseq_from_seq1.fasta"
                refseqList = []
                for marker in markerList:
                    fastaFile = f"{args.output}/{marker}/{marker}.uniqueSeqs.fasta"
                    firstSeq = list(SeqIO.parse(fastaFile, "fasta"))[0]
                    firstSeq.id = marker
                    firstSeq.name = marker
                    firstSeq.description = ""
                    refseqList.append(firstSeq)

                wsh = open(refseqFile, "w")
                SeqIO.write(refseqList, wsh, "fasta")
                wsh.close()

            match = 2
            mismatch = -1
            scoring = swalign.NucleotideScoringMatrix(match, mismatch)
            sw = swalign.LocalAlignment(scoring)

            refseqExist = set()
            markerDirList = []

            for record in SeqIO.parse(refseqFile, "fasta"):
                markerName = record.id
                refseqExist.add(markerName)
                seqStr = str(record.seq)
                seqLen = len(seqStr)
                refseq_5p25nt = seqStr[0:25]
                refseq_3p25nt = seqStr[seqLen-25:seqLen]
                markerDirList.append((markerName,refseq_5p25nt, refseq_3p25nt, args.filter, args.alnpct))
  
            for m in markerList:
                if m not in refseqExist:
                    print(f"Error: Reference sequence fasta file {refseqFile} does not contain reference sequence for marker {m}!")
                    sys.exit()

            pool = multiprocessing.Pool(processes= args.job)
            pool.starmap(finalProcess, markerDirList)
            pool.close()

    if ("4" not in args.skip):
        logging.info("Step 4: delete intermediate")
        cmd = f"rm {args.output}/*/*fastq.gz"
        returned_value = subprocess.call(cmd, shell=True)
        cmd = f"rm {args.output}/*/*.RData"
        returned_value = subprocess.call(cmd, shell=True)



def splitByCutadapt(sampleName, file1, file2):
    try:
        sampleDir = f"{args.output}/{sampleName}"
        if (not os.path.exists(sampleDir)):
            os.mkdir(sampleDir)

        #run cutadapt to demultiplexing by primers
        cmd = f"{cutadaptCMD} --quiet -e {args.primerErrorRate} --minimum-length={args.minHaplotypeLength} --discard-untrimmed -g file:{primerFileF} -G file:{primerFileR} -o {sampleDir}/{{name1}}_{{name2}}_R1.fastq -p {sampleDir}/{{name1}}_{{name2}}_R2.fastq {file1} {file2} "
        returned_value = subprocess.call(cmd, shell=True)
        logging.info(f"the command is {cmd}")
        logging.info(f"demultiplexing {sampleName} done: {returned_value}")


        # move splitted file to corresponding marker directory
        for markerName in markerList:
            if (os.path.exists(f"{sampleDir}/{markerName}_{markerName}_R1.fastq")):
                cmd1 = f"mv {sampleDir}/{markerName}_{markerName}_R1.fastq {args.output}/{markerName}/{sampleName}_R1.fastq"
                cmd2 = f"mv {sampleDir}/{markerName}_{markerName}_R2.fastq {args.output}/{markerName}/{sampleName}_R2.fastq"

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
    returned_value =-99
    try:
        cmd = f"{dadaRscript} {args.output}/{markerPath}"
        returned_value = subprocess.call(cmd, shell=True)
    except Exception as e:
        print(f'Caught exception in job {markerPath} ')
        traceback.print_exc()
        print()
        raise e
    return returned_value


def finalProcess(markerPath, refseq_5p25nt, refseq_3p25nt, alignmethod, alignpct):

    try:
        BLAST5prime = set()
        BLAST3prime = set()
        seq2ID = {}
        acceptedIds = []
        acceptedRecords= []
        removedRecords= []

        fastaFile = f"{args.output}/{markerPath}/{markerPath}.uniqueSeqs.fasta"
        fastaFilteredFile = f"{args.output}/{markerPath}/{markerPath}.uniqueSeqs.filtered.fasta"
        fastaRemovedFile = f"{args.output}/{markerPath}/{markerPath}.uniqueSeqs.removed.fasta"
        blastResult = f"{args.output}/{markerPath}/blastresults"

        if (os.path.isfile(fastaFile)):
            if (alignmethod == "2"):
                matchedNT = alignpct*25 
                for record in SeqIO.parse(fastaFile, "fasta"):
                    seq2ID[str(record.seq)] = record.id
                    markerName = record.id
                    seqStr = str(record.seq)
                    seqLen = len(seqStr)
                    query_5p25nt = seqStr[0:25]
                    query_3p25nt = seqStr[seqLen-25:seqLen]
                    alignment1 = sw.align(refseq_5p25nt, query_5p25nt)
                    alignment2 = sw.align(refseq_3p25nt, query_3p25nt)
                    if (alignment1.matches>matchedNT and alignment2.matches>matchedNT):
                        acceptedIds.append(record.id)
                        acceptedRecords.append(record)
                    else:
                        removedRecords.append(record)        
                wsh = open(fastaFilteredFile, "w")
                wsh2 = open(fastaRemovedFile, "w")
                SeqIO.write(acceptedRecords, wsh, "fasta")
                SeqIO.write(removedRecords, wsh2, "fasta")
                wsh.close()
                wsh2.close()
            if (alignmethod == "3"):
                ref5Query = f"{args.output}/{markerPath}/{markerPath}.ref5.fasta"
                ref3Query = f"{args.output}/{markerPath}/{markerPath}.ref3.fasta"
                blastOut5 = f"{args.output}/{markerPath}/{markerPath}.ref5.blast"
                blastOut3 = f"{args.output}/{markerPath}/{markerPath}.ref3.blast"
                w1 = open(ref5Query, "w")
                w1.write(">ref5\n$refseq_5p25nt\n")
                w1.close()
                w1 = open(ref3Query, "w")
                w1.write(">ref5\n$refseq_3p25nt\n")
                w1.close()
                cmd = f"{makeblastdbCmd} -in {fastaFile} -dbtype nucl"
                subprocess.call(cmd, shell=True)
                cmd = f"{blastCmd} -query {ref5Query} -db {fastaFile} -task \"blastn-short\" -outfmt 6 -out {blastOut5}"
                subprocess.call(cmd, shell=True)
                cmd = f"{blastCmd} -query {ref3Query} -db {fastaFile} -task \"blastn-short\" -outfmt 6 -out {blastOut3}"
                subprocess.call(cmd, shell=True)

                hits5 = []
                hits3 = []
                with open (blastOut5, "r") as BLASTIN:
                    for line in BLASTIN:
                        dataFields = line.split("\t")
                        if (float(dataFields[2])/100 >= alignpct ):
                            hits5.append(dataFields[1])
                with open (blastOut3, "r") as BLASTIN:
                    for line in BLASTIN:
                        dataFields = line.split("\t")
                        if (float(dataFields[2])/100 >= alignpct ):
                            hits3.append(dataFields[1])
                acceptedIds = [value for value in hits5 if value in hits3]

                for record in SeqIO.parse(fastaFile, "fasta"):
                    if(record.id in acceptedIds):
                        acceptedRecords.append(record)
                    else:
                        removedRecords.append(record)                   
                wsh = open(fastaFilteredFile, "w")
                wsh2 = open(fastaRemovedFile, "w")
                SeqIO.write(acceptedRecords, wsh, "fasta")
                SeqIO.write(removedRecords, wsh2, "fasta")
                wsh.close()
                wsh2.close()

        ## filter .seqtab.csv file
        csvFile = f"{args.output}/{markerPath}/{markerPath}.seqtab.csv"
        modCsvFile = f"{args.output}/{markerPath}/{markerPath}.seqtab.mod.csv"
        oriCsvFile = f"{args.output}/{markerPath}/{markerPath}.seqtab.ori.csv"

        dataMatrix= pd.read_csv(csvFile, header=0)
        samplesInFile = [item.replace("_F_filt.fastq.gz","") for item in dataMatrix.columns]
        samplesInFile[0] = "alleleSequence"
        dataMatrix.columns = samplesInFile
        dataMatrix['Id']= dataMatrix["alleleSequence"].map(seq2ID)

        newColumnName = sampleList
        newColumnName[:0] = ["Id", "alleleSequence"]

        oriMatrix_redexed = dataMatrix.reindex(columns=newColumnName, fill_value=0)  
        oriMatrix_redexed.to_csv(oriCsvFile, quotechar='"', quoting=csv.QUOTE_NONNUMERIC, index=False)

        filteredDataMatrix = oriMatrix_redexed.loc[dataMatrix['Id'].isin(acceptedIds)]
        filteredDataMatrix.to_csv(modCsvFile, quotechar='"', quoting=csv.QUOTE_NONNUMERIC, index=False)    
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


