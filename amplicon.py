#!/usr/bin/python3

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


from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Alphabet import IUPAC, Gapped
from Bio import SeqIO

#import uuid
#from scipy import stats
#import numpy as np
#import itertools



## dependency on BioHPC
#export PATH=/programs/muscle:$PATH:/programs/bbmap-38.45:$PATH

#some constants
haplotypeAvgToLengthRatio =2 #For each marker, if an allele has less than 1/haplotypeAvgToLengthRatio of the averagelength, this allele is skipped

def main():

    commandLine = " ".join(sys.argv[1:])

    global args
    global sampleList
    global sampleToFileList
    global markerList
    global primerFile
    global tagBySampleDir
    global markerDir
    global HapGTFile
    global HapMatrixFile
    global modDir
    global correctErrorFileList
    global readCountMatrixFile
    global alleleFastaFile
    global topAlleleFastaFile

    #check dependencies
    print ("Checking dependencies:")
    if (not checkApp("bbmerge.sh")):
        sys.exit()

    if (not checkApp("cutadapt")):
        sys.exit()

    if (not checkApp("muscle")):
        sys.exit()


    parser = argparse.ArgumentParser(description='Run GATK Haplotype Caller on ampseq data.')

    # Required arguments
    parser.add_argument('-s','--sample',type=str,required=True,help='Sample file. Tab delimited text file with 3 or 4 columns: sample_Name, plate_well, fastq_file1, (optional)fastq_file2. Plate_well is a string to uniquely define sample if sample names are duplicated.')
    parser.add_argument('-k','--key',type=str,required=True,help='Key file. Tab delimited text file with 3 columns: amplicon name, 5 prime sequence, 3 prime sequence')
    parser.add_argument('-o','--output',type=str,required=True,help='Output directory')


    # Optional arguments
    parser.add_argument('-i','--skip',type=str,required=False,default="",help='Skip steps. e.g. "-i 12" to skip steps 1 and 2. the steps are: 1. split reads by primers; 2. identify haplotypes across population, and optionally run PCR error correction if set "-e 1"; 3 call genotypes')
    parser.add_argument('-j','--job',type=int,required=False,default=8,help='Number of simultaneous jobs. Default:8')
    parser.add_argument('-t','--thread',type=int,required=False,default=1,help='Number of threads per job. Default:1')
    parser.add_argument('-c','--minSamplePerHaplotype',type=int,required=False,default=10,help='Minimum number of samples per haplotypes. Default:10')
    parser.add_argument('-n','--maxHaplotypePerSample',type=int,required=False,default=20,help='Maximum number of unique haplotypes per sample in the first pass, no matter what the ploidy level of the individual. Default:20')
    parser.add_argument('-m','--maxHaplotypeInPopulation',type=int,required=False,default=1000,help='Maximum number of haplotypes per marker in the population. Default:1000')
    parser.add_argument('-a','--maf',type=float,required=False,default=0.01,help='Minimum minor allele frequency Default:0.01')
    parser.add_argument('-l','--minHaplotypeLength',type=int,required=False,default=50,help='Minimum haplotype length (after removing the primers. It must be an integer 1 or larger.) Default:50')
    parser.add_argument('-d','--mergeDuplicate',type=int,required=False,default=1,help='Whether to merge the duplicate samples. 1: merge; 0: not merge and the duplicated sample will be named <sampleName>__<index starting from 1> . Default:1')
    parser.add_argument('-e','--PCRErrorCorr',type=int,required=False,default=0,help='Correct PCR errors based on allele frequency (only applicable for biparental families). 0: No correction; 1: Correct error in bi-parental population based on allele read count distribution in the population. Default:0, no correction')
    parser.add_argument('-p','--ploidy',type=int,required=False,default=2,help='Ploidy, default 2')
    parser.add_argument('-r','--maxAlleleReadCountRatio',type=int,required=False,default=20,help='Maximum read count ratio between the two alleles in each sample, default 20')


    if sys.version_info[0] < 3:
        raise Exception("This code requires Python 3.")

    args=parser.parse_args()


    sampleList = []
    sampleToFileList = []
    markerList= []
    correctErrorFileList = []

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
    markerDir = f"{args.output}/markerDir"
    HapGTFile = f"{args.output}/hap_genotype"
    HapMatrixFile = f"{args.output}/hap_genotype_matrix"
    readCountMatrixFile = f"{args.output}/markerToSampleReadCountMatrix"
    alleleFastaFile = f"{args.output}/HaplotypeAllele.fasta"
    topAlleleFastaFile = f"{args.output}/topHaplotypeAllele.fasta"
    modDir = f"{args.output}/mod1"   # store data from error correction 1


    if (not os.path.exists(args.output)):
        os.mkdir(args.output)
    if (not os.path.exists(tagBySampleDir)):
        os.mkdir(tagBySampleDir)
    if (not os.path.exists(markerDir)):
        os.mkdir(markerDir)
    if ((args.PCRErrorCorr ==1) and (not os.path.exists(modDir))):
        os.mkdir(modDir)


    logging.basicConfig(filename=args.output + '/run.log',level=logging.DEBUG)

    logging.debug(f"Command: {commandLine}")

    # Process the key file, and prepare linked adapter as required by cutadapt
    primerfh = open(primerFile, "w")
    with open(args.key, 'r') as fhk:
        for line in fhk:
            if (not re.match("\w", line)):
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
        fhk.close()
    primerfh.close()

    #process sample file, currently, only paired end reads are supported
    ## first check and merge duplicate samples
    checkSampleDup = {}
    fileMerged = {}

    with open(args.sample, 'r') as fhs:
        for line in fhs:
            if (not re.match("\w", line)):
                continue
            line = line.rstrip()
            fieldArray = line.split(sep="\t")
            sampleName = re.sub("\s", "", fieldArray[0])
            plateWell = re.sub("\s", "", fieldArray[1])
            if (sampleName in checkSampleDup):
                checkSampleDup[sampleName].append((fieldArray[2], fieldArray[3]))
            else:
                checkSampleDup[sampleName] = [(fieldArray[2], fieldArray[3])]
        fhs.close()

    # execute file merging only if step 1 not skipped

    for sampleName, files in checkSampleDup.items():
        if (len(files)>1):
            if (args.mergeDuplicate == 1):
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
            else:
                fileMerged[sampleName] = "duplicated"


    with open(args.sample, 'r') as fhs:
        for line in fhs:
            if (not re.match("\w", line)):
                continue
            line = line.rstrip()
            fieldArray = line.split(sep="\t")
            if (len(fieldArray) < 4):
                print(f"Error: Sample file must have at least four columns, and with no header line. Single-end reads are not supported now. Will be added later")
                sys.exit()

            sampleName = re.sub("\s", "", fieldArray[0])
            plateWell = re.sub("\s", "", fieldArray[1])
            if (sampleName in fileMerged):
                if ((args.mergeDuplicate == 1) and (fileMerged[sampleName] == "done")):
                    continue
                elif (args.mergeDuplicate == 1):
                    fieldArray[2] = fileMerged[sampleName][0]
                    fieldArray[3] = fileMerged[sampleName][1]
                    fileMerged[sampleName] = "done"
                else:
                    sampleName = sampleName + "$" + plateWell

            sampleList.append(sampleName)

            if (not os.path.isfile(fieldArray[2])):
                print(f"Error: Sample fastq file {fieldArray[2]} does not exist!")
                sys.exit()
            if (not os.path.isfile(fieldArray[3])):
                print(f"Error: Sample fastq file {fieldArray[3]} does not exist!")
                sys.exit()
                
            sampleToFileList.append((sampleName, fieldArray[2], fieldArray[3]))

    # split the read by primers, and remove primer from each read, and collapase identical reads, and keep top <arg.maxHaplotypePerSample> tags per sample
    if ("1" not in args.skip):
        logging.info("Step 1: splitByPrimer")
        splitByPrimer()

    # get sequence tag list across population
    if ("2" not in args.skip):
        logging.info("Step 2: getTagList")
        getTagList()
        if (args.PCRErrorCorr == 1):
            correctPCRError()

    # callHapGenotypes
    if ("3" not in args.skip):
        logging.info("Step 3: call haplotype genotypes")
        callHapGenotypes()

    # convert to LepMap output
    #logging.info("Step 4: output LepMap file")
    #toLepMap()


def splitByPrimer():
    print("Run splitByPrimer ...")
    pool = multiprocessing.Pool(processes= args.job)
    pool.starmap(splitByCutadapt, sampleToFileList)
    pool.close()

    mmFh = open (readCountMatrixFile, "w")
    mmFh.write("\t")
    mmFh.write("\t".join(markerList))
    mmFh.write("\n")
    mmFh.close()
    os.system(f"cat {tagBySampleDir}/*.readcount >> {readCountMatrixFile}")
    os.system(f"rm {tagBySampleDir}/*.readcount")




def splitByCutadapt(sampleName, file1, file2):
    try:
        sampleDir = f"{args.output}/{sampleName}"
        if (not os.path.exists(sampleDir)):
            os.mkdir(sampleDir)

        #contig the paired end reads
        cmd = f"bbmerge.sh t={args.thread} in1={file1} in2={file2} outm={sampleDir}/contig.fastq"
        logging.info(f"Process {sampleName}: {cmd}")
        returned_value = subprocess.call(cmd, shell=True)
        logging.info(f"contiging {sampleName} done: {returned_value}")

        #run cutadapt to demultiplexing by primers
        cmd = f"cutadapt --quiet -e 0.1 --minimum-length={args.minHaplotypeLength} --trimmed-only -g file:{primerFile} -o {sampleDir}/{{name}}.fastq {sampleDir}/contig.fastq "
        returned_value = subprocess.call(cmd, shell=True)
        logging.info(f"demultiplexing {sampleName} done: {returned_value}")

        #collapse identical reads
        tagBySampleFile = f"{tagBySampleDir}/{sampleName}.tbs"
        tbsfh = open(tagBySampleFile, "w")

        readCountMatrix = {}
        
        rcFh = open (f"{tagBySampleDir}/{sampleName}.readcount", "w")
        rcFh.write(f"{sampleName}")
        for marker in markerList:
            readCount =0
            markerFile = f"{sampleDir}/{marker}.fastq"
            collapsedMarkerFile = f"{sampleDir}/{marker}.collapsed"
            if (os.path.exists(markerFile)):
                cmd=f"awk 'NR%4==2' {markerFile}|LC_ALL=C sort |uniq -c |LC_ALL=C  sort -k1,1rn > {collapsedMarkerFile}"
                returned_value = subprocess.call(cmd, shell=True)
                tagCount = 0
                topAlleleReadCount = 0
                
                with open(collapsedMarkerFile, 'r') as cfh:
                    for line in cfh:
                        [copyNumber, seqStr] = line.split()
                        copyNumber = int(copyNumber)
                        readCount += copyNumber
                        if (tagCount==0):
                            topAlleleReadCount = copyNumber
                        if (tagCount < args.maxHaplotypePerSample):
                            if ((topAlleleReadCount/copyNumber) < args.maxAlleleReadCountRatio):
                                tagCount+=1
                                tbsfh.write(f"{sampleName}\t{marker}\t{seqStr}\t{copyNumber}\n")
                    cfh.close()
            rcFh.write(f"\t{readCount}")
        tbsfh.close()
        rcFh.write("\n")
        rcFh.close()
        logging.info(f"collapsing contig {sampleName} done")
        os.system("rm -rf %s" % sampleDir)
        return 1 
    except Exception as e:
        print(f'Caught exception in job {sampleName} {marker}')
        traceback.print_exc()
        print()
        raise e

def getTagList():
    print ("Collapsing tags across samples")
    
    #collapse tags across samples
    # get sample and treat count for each unique sequence tag
    cmd = f"cut -f 2,3,4  {tagBySampleDir}/*.tbs"
    proc = subprocess.Popen(cmd,stdout=subprocess.PIPE,shell=True)
    tagToReadCount = {}
    markerTotalLen = {}
    markerSeqStrcount = {}

    for marker in markerList:
        markerTotalLen[marker]=0
        markerSeqStrcount[marker]=0

    with proc.stdout as f:
        for line in f:
            [marker, seqStr, copyNumber]=line.split()

            marker = marker.decode('utf-8')
            seqStr = seqStr.decode('utf-8')
            copyNumber = int(copyNumber.decode('utf-8'))

            markerTotalLen[marker] += len(seqStr)
            markerSeqStrcount[marker] +=1

            if ('N' in seqStr):
                continue

            if (marker not in tagToReadCount):
                tagToReadCount[marker] = {}
            if (seqStr not in tagToReadCount[marker]):
                tagToReadCount[marker][seqStr] = [1, copyNumber]
            else:
                tagToReadCount[marker][seqStr][0] += 1
                tagToReadCount[marker][seqStr][1] += copyNumber
        f.close()
    #write the results to one file per marker, reverse sorted by sample count per sequence tag
    for marker in markerList:
        markerAvgLen = 0
        if (markerSeqStrcount[marker]>0):
            markerMinLen = int(markerTotalLen[marker]/markerSeqStrcount[marker]/haplotypeAvgToLengthRatio)

        if (marker in tagToReadCount):
            tagId=0
            firstPassFile = f"{markerDir}/{marker}.firstpass"
            modOutFile = f"{modDir}/{marker}"
            tagTableFh = open(firstPassFile, "w")
            correctErrorFileList.append((firstPassFile, modOutFile))
            #duelCount: each seqstr has two counts stored as a list: by r by sample and by reads
            for (seqStr, duelCount) in sorted(tagToReadCount[marker].items(), key=lambda x: x[1][0], reverse=True):
                if (len(seqStr) < markerMinLen):
                    continue
                if ((duelCount[0] >=args.minSamplePerHaplotype) and (tagId < args.maxHaplotypeInPopulation)):
                    tagId +=1
                    tagTableFh.write(f"{tagId}\t{seqStr}\t{duelCount[0]}\t{duelCount[1]}\n");
            tagTableFh.close()

    tagToReadCount.clear()
    

def correctPCRError():
    print("Run PCR error correction ...")
    pool = multiprocessing.Pool(processes= args.job)
    pool.starmap(main_collapse, correctErrorFileList)
    pool.close()
    os.system(f"rm {modDir}/*.fa")
    #os.system(f"rm {modDir}/*.aln")

def callHapGenotypes():
    print ("Call Genotypes")
    #collect all valid haplotype alleles with their ID
    seqTagToId = {}
    GTMatrix = {}
    for marker in markerList:
        GTMatrix[marker] = {}
        if (args.PCRErrorCorr == 1):
            #no use error correction method 1 (-e 1)
            if (os.path.isfile(f"{modDir}/{marker}.mod")):
                with open(f"{modDir}/{marker}.mod", 'r') as ms:
                    for line in ms:
                        [oriTagId, seqStr, xx, yy, informaticeSites, modTagId]  = line.split()
                        if marker not in seqTagToId:
                            seqTagToId[marker] = {}
                        seqTagToId[marker][seqStr] = int(modTagId)
                    ms.close()
        else:
            #use firstpass (-e 0)
            if (os.path.isfile(f"{markerDir}/{marker}.firstpass")):
                with open(f"{markerDir}/{marker}.firstpass", 'r') as ms:
                    for line in ms:
                        [tagId, seqStr, sampleNumber, copyNumber]  = line.split()
                        if marker not in seqTagToId:
                            seqTagToId[marker] = {}
                        seqTagToId[marker][seqStr] = int(tagId)
                    ms.close()

    #get genotypes matrix
    for sampl in sampleList:
        if (os.path.isfile(f"{tagBySampleDir}/{sampl}.tbs")):
            markerTagIDToReadCount = {}
            with open(f"{tagBySampleDir}/{sampl}.tbs", 'r') as sfh:
                for line in sfh:
                    [sss, marker, seqStr, copyNumber]  = line.split()
                    copyNumber = int(copyNumber)
                    if ((marker in seqTagToId) and (seqStr in seqTagToId[marker])):
                        alleleId = seqTagToId[marker][seqStr]
                        if (marker not in markerTagIDToReadCount):
                            markerTagIDToReadCount[marker] = {}
                        if (alleleId in markerTagIDToReadCount[marker]):
                            markerTagIDToReadCount[marker][alleleId] += copyNumber
                        else:
                            markerTagIDToReadCount[marker][alleleId] = copyNumber
                sfh.close()
                ## within the sample, get top two allele per marker
                for marker in markerTagIDToReadCount:
                    sorted_alleles = sorted(markerTagIDToReadCount[marker].items(), key=lambda x: x[1], reverse=True)
                    alleleCounter=0
                    for (alleleid, copies) in sorted_alleles:
                        if (alleleCounter < args.ploidy):
                            alleleCounter +=1
                            if (sampl in GTMatrix[marker]):
                                GTMatrix[marker][sampl].append((alleleid, copies))
                            else:
                                GTMatrix[marker][sampl] = [(alleleid, copies) ]

    #output matrix
    gtfh = open(HapGTFile, "w")
    gtfh.write ("Locus\tHaplotypes\t")
    gtfh.write ("\t".join(sampleList))
    gtfh.write ("\n")


    outputMarkerList = []
    outputGT = {}
    markerToFinalAllele = {}
    for marker in markerList:
        if (marker not in GTMatrix):
            continue
        outputGT[marker] = {}
        markerToFinalAllele[marker] = []
        #gather allele frq, and only keep allele with maf > args.maf
        alleleFrq = {} ## store allele frequency for 
        total = 0   ## total number of genotyped gametes in the population for this marker, for calculate allele frequency
        keptAlleles = {}  ## next block of code calculate allele frq for each allele, only allele above threshold are in keptallelels
        for sampl in sampleList:
            if (sampl in GTMatrix[marker]):
                alleles = GTMatrix[marker][sampl]
                total += args.ploidy
                if (len(alleles) == 1):
                    (alleleid, copyNumber) = alleles[0]
                    if (alleleid in alleleFrq):
                        alleleFrq[alleleid] +=args.ploidy
                    else:
                        alleleFrq[alleleid] = args.ploidy               
                else:
                    alleleCounter=0
                    for (alleleid, copies) in alleles:
                        alleleCounter+=1
                        if (alleleid in alleleFrq):
                            alleleFrq[alleleid] +=1
                        else:
                            alleleFrq[alleleid] =1


        if (total==0):
            continue


        for aid, copyNumber in alleleFrq.items():
            frq = copyNumber/total
            if (frq > args.maf):
                keptAlleles[aid] = 1

        if (len(keptAlleles) == 0):
            continue

        # write marker output, only use keptAlleles
        gtline = ""
        alleleFrq = {}
        total = 0 

        #gather allele frq
        for sampl in sampleList:
            if (sampl in GTMatrix[marker]):             
                alleles = GTMatrix[marker][sampl]
                total += args.ploidy
                alleleCounter=0
                gtline_p1 = []
                gtline_p2 = []
                for (alleleid, copies) in alleles:
                    if (alleleid in keptAlleles):
                        alleleCounter += 1
                        gtline_p1.append(str(alleleid))
                        gtline_p2.append(str(copies))

                if (len(gtline_p1) == 0):
                    gtline+="\t./.:0"
                elif (len(gtline_p1) == 1):
                    if (gtline_p1[0] in alleleFrq):
                        alleleFrq[gtline_p1[0]] +=args.ploidy
                    else:
                        alleleFrq[gtline_p1[0]] =args.ploidy
                    if (args.ploidy==2):
                        gtline+= f"\t{gtline_p1[0]}/{gtline_p1[0]}:{gtline_p2[0]}"
                        outputGT[marker][sampl] = f"{gtline_p1[0]}|{gtline_p1[0]}"
                    else:
                        gtline+=f"\t{gtline_p1[0]}:{gtline_p2[0]}"
                        outputGT[marker][sampl] = f"{gtline_p1[0]}"
                else:
                    for alleleid in gtline_p1:
                        if (alleleid in alleleFrq):
                            alleleFrq[alleleid] +=1
                        else:
                            alleleFrq[alleleid] =1
                    p1 =  "/".join(gtline_p1)
                    p2 = ",".join(gtline_p2)
                    gtline+= f"\t{p1}:{p2}"
                    outputGT[marker][sampl] = p1.replace("/", "|")
            else:
                gtline += "\t./.:0"

        frqline = ""
        if (total==0):
            continue
        alleleFrqSorted = sorted(alleleFrq.items(), key=lambda x: x[1], reverse=True)
        for aaa in alleleFrqSorted:
            (aid, copyNumber) = aaa
            frq = copyNumber/total
            frqline += f"{aid}({frq:4.2f});"
            markerToFinalAllele[marker].append(aid)
        gtfh.write (f"{marker}\t{frqline}{gtline}\n")
        outputMarkerList.append(marker)
    gtfh.close()

    ## write haplotype to sample matrix
    hmFh = open (HapMatrixFile, "w")
    hmFh.write("\t")
    hmFh.write("\t".join(outputMarkerList))
    hmFh.write("\n")
    for sampl in sampleList:
        hmFh.write(sampl)
        for marker in outputMarkerList:
            if ((marker in outputGT) and (sampl in outputGT[marker] )):
                hmFh.write(f"\t{outputGT[marker][sampl]}")
            else:
                hmFh.write("\t.|.")
        hmFh.write("\n")
    hmFh.close()


    ## write a fasta file with all alleles, and a fasta file with only the top allele
    fastaFh = open (alleleFastaFile, "w")
    topFh = open (topAlleleFastaFile, "w")
    for marker in outputMarkerList:
        if marker in markerToFinalAllele:
            alleleToSeq = {}
            with open(f"{markerDir}/{marker}.firstpass", 'r') as ms:
                for line in ms:
                    [tagId, seqStr, sampleNumber, copyNumber]  = line.split()
                    alleleToSeq[tagId] =seqStr
            topAllele= markerToFinalAllele[marker][0]
            topFh.write(f">{marker}\n{alleleToSeq[topAllele]}\n")
            for allele in markerToFinalAllele[marker]:
                fastaFh.write(f">{marker}#{allele}\n{alleleToSeq[allele]}\n")

    fastaFh.close()


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


########### the following code is for PCR error correction ##############################################

alphabet = Gapped(IUPAC.ambiguous_dna)

def align_muscle(inputfa):

	from Bio.Align.Applications import MuscleCommandline	
	muscle_cline = MuscleCommandline(input=inputfa)
	muscle_cline.set_parameter("quiet", True)
	stdout, stderr = muscle_cline()
	from io import StringIO
	from Bio import AlignIO
	
	id_order=[]
	for record in SeqIO.parse(inputfa, "fasta"):
    		id_order.append(record.id)

	align_dict = SeqIO.to_dict(SeqIO.parse(StringIO(stdout), "fasta"))

	sort_align = [] # Setup an empty list
	for record in id_order:
		sort_align.append(align_dict[record])
	SeqIO.write(sort_align, inputfa+".aln", "fasta")
	align = AlignIO.read(inputfa+".aln", "fasta")
	#print(align)
	return align

def sliding_entropy(alignment):
	summary_align = AlignInfo.SummaryInfo(alignment)
	
	consensus = summary_align.gap_consensus(threshold=0.01)
	#print(consensus)
	entropy_list=[]
	for i in range(len(consensus)-5):
		p=cal_entropy(consensus[i:i+5])
		#print(p)
		entropy_list.append(p)
	for i in range(len(consensus)-5, len(consensus)):
		p=cal_entropy(consensus[i:])
		#print(p)
		entropy_list.append(p)

	return entropy_list

def cal_entropy(test_str):
	### get entropy for each position

	res = dict(Counter(test_str))
	
	#print(res)
	px=[x/float(sum(list( res.values()))) for x in list(res.values())]
	echo_entropy=[ x*log(x) for x in px]
	Shannon_column_entropy = - sum(echo_entropy)
	return Shannon_column_entropy



def cal_MAF(align, freq_list):
	### get MAF for each position
	maf_list=[]
	#print(align.format("phylip"))
	align_len=align.get_alignment_length()
	for i in range(align_len):
		#print(len(align),i)
		test_str=align[:, i]
		#print(align[:, i])
		res = Counter(test_str)
		
		s = [(k, res[k]) for k in sorted(res, key=res.get, reverse=True)]
		major=s[0][0]
		#print(major)
		if 1:
			occur_list= [  1 if test_str[x]== major else 0  for x in range(len(align) )]
			#print(occur_list)
			multi_count= [  freq_list[i]*occur_list[i]   for i in range(len(align))]
			#print(multi_count)
			#print(sum(freq_list))
			f=sum(multi_count)/float(sum(freq_list))
			maf=1-f
		
		maf_list.append(maf)
	#print(maf_list)
	return maf_list

def len_fil(seqdic,outfile):
	###############################	
	##  < 50% is not allowed. 
	##
	keylist=list(seqdic.keys())
	
	len_dic=dict(zip(keylist,[len(seqdic[i][0]) for i in keylist ]))
	#print(len_dic)
	maxLen=max(list(len_dic.values()))
	occ_list=[]
	
	oup=open(outfile,"w")
	kept_index=[]
	for i in keylist:
		if len_dic[i]/float(maxLen)>0.5:
			oup.write(">%s\n%s\n" % (i,seqdic[i][0]) )
			occ_list.append(seqdic[i][1] )
			kept_index.append(i)
	oup.close()
	#print(occ_list)
	return [occ_list,set(kept_index)]	

def main_collapse(inp,ouf):

	print("%s\n" % "#########################################################")
	print(inp)
	##parse the input file and get fasta and counts for each seq.
	frq=[]
	#tmp=str(uuid.uuid4())+".fa"
	tmp=ouf+".fa"

	##get haplotype from file

	file = Path(inp)
	size = file.stat().st_size
	if size<2:
		return "firstpass empty file"

	
	seq_dic={}
	for line in open(inp,"r"):
		L=line.strip().split()
		seq_dic[L[0]]=[L[1],int(L[3])]

	##########################
	### length filtering
	### Here I used two cirteria
	###  1.abs(zscore) <1
	###  2. frq of this haplotype <0.05

	[frq,kept_index]=len_fil(seq_dic,tmp)	
	#print(frq,kept_index)
	file = Path(tmp)
	size = file.stat().st_size
	if size<2:
		return "aftering length filtering empty file"
		

	align=align_muscle(tmp)
	maf=cal_MAF(align, frq)
	entropy=sliding_entropy(align)
	#print(maf)
	#print(entropy)
	#print(', '.join(f'{x:.3f}' for x in maf))
	#print(', '.join(f'{x:.3f}' for x in entropy))
	good_site=[]
	for i in range(len(maf)):
		if abs(maf[i]-0.25)<0.10:
			#print(maf[i])
			start=max(0,i-5)
			end=min(i+2,len(maf) )
			print(start,end)
			print(entropy[start:end])
			minentropy=min(entropy[start:end],default=entropy[0])
			if minentropy>0.3:
				good_site.append(i)
		if abs(maf[i]-0.50)<0.10:
			#print(maf[i])
			start=max(0, i-5)
			end=min(i+2,len(maf) )
			print(start,end)
			print(entropy[start:end])
			minentropy=min(entropy[start:end],default=entropy[0])
			if minentropy>0.3:
				good_site.append(i)
	#print(good_site)

	if len(good_site)==0:
		from operator import itemgetter
		indices, L_sorted = zip(*sorted(enumerate(maf), key=itemgetter(1), reverse=True))
		for i in indices:
			if i==0:
				i=1
			#print(i)
			start=max(0, i-5)
			end=min(i+2,len(maf) )
			minentropy=min(entropy[start:end],default=entropy[0])
			if minentropy>0.3:
				good_site.append(i)
				break

	if len(good_site)==0:
		from operator import itemgetter
		indices, L_sorted = zip(*sorted(enumerate(maf), key=itemgetter(1), reverse=False))
		good_site.append(indices[0])

	for i in good_site:
		print(align[:, i])

	oup=open(ouf+".mod","w")
	indexout_haplotype={}
	n=1
	#print(inp)

	align_dic={}
	for record in align:
		align_dic[record.id]=record.seq
	print(inp)
	with open(inp) as myFile:
		for num, line in enumerate(myFile, 0):
			L=line.strip().split("\t")
			#print(L[0])
			if L[0] in kept_index:
				print(L[0])
				out_str= "".join([align_dic[L[0]][i]  for i in good_site])
				if out_str not in indexout_haplotype:
					indexout_haplotype[out_str]=n
					n+=1
				
				oup.write("%s\t%s\t%s\n"  % ( line.strip(), out_str, indexout_haplotype[out_str])   )
	oup.close()
	#os.system("rm %s" % tmp)
	#os.system("rm %s.aln" % tmp )



########### end of  PCR error correction code #################################

if __name__=="__main__":
    main()


