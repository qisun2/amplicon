#!/usr/bin/env python3

import logging
from os import listdir
from os.path import isfile, join
import argparse
import pathlib
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


#some constants
bwa = "bwa"
bbmerge = "/programs/bbmap-38.45/bbmerge.sh"
gatk = "/programs/gatk4/gatk"
samtools = "samtools"
picard = "java -jar /programs/picard-tools-2.19.2/picard.jar"
RELEASE = "/programs/sentieon-genomics-201808.05"
os.environ["SENTIEON_LICENSE"] = "cbsulogin2.tc.cornell.edu:8990"
slurmScript = "/home/vitisgen/tools/run_gatk.py"

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
    global slurmSampleList
    global slurmBatchFile


    parser = argparse.ArgumentParser(description='Run GATK Haplotype Caller on ampseq data, using slurm.')

    # Required arguments
    parser.add_argument('-s','--sample',type=str,required=True,help='Sample file. Tab delimited text file with 3 or 4 columns: sample_Name, plate_well, fastq_file1, (optional)fastq_file2. Plate_well is a string to uniquely define sample if sample names are duplicated.')
    parser.add_argument('-o','--output',type=str,required=True,help='Output directory')
    parser.add_argument('-r','--reference',type=str,required=True,help='directory of reference db. the directory should be prepared with the script /home/vitisgen/tools/gatk_genome.sh ref.fasta genomeDir')
    # Optional arguments
    parser.add_argument('-m','--trim',type=int,required=False,default="25",help='Trim both the 5- and 3-prime end of the amplicon by this length. Default 25. Set to 0 for no trimming. The is necessary as primer regions are not reliable for genotyping.')
    parser.add_argument('-j','--job',type=int,required=False,default=8,help='Number of simultaneous jobs. Default:8')
    parser.add_argument('-t','--thread',type=int,required=False,default=1,help='Number of threads per job. Default:1')
    parser.add_argument('-d','--mergeDuplicate',type=int,required=False,default=1,help='Whether to merge the duplicate samples. 1: merge; 0: not merge and the duplicated sample will be named <sampleName>__<index starting from 1> . Default:1')
    parser.add_argument('-x','--slurmcluster',type=str,required=False,default="",help='Slurm cluster name. Slurm is configured to run on Cornell BioHPC. e.g. "-x cbsumm10"')
    parser.add_argument('-y','--slurmBatchSize',type=int,required=False,default=10,help='Slurm job size. Number of samples per slurm job. default 10')
    parser.add_argument('-i','--skip',type=str,required=False,default="",help='Skip steps. e.g. "-i 1" to skip steps 1. the steps are: 1. haplotype caller; 2. genotype')


    if sys.version_info[0] < 3:
        raise Exception("This code requires Python 3.")

    args=parser.parse_args()


    sampleList = []
    sampleToFileList = []
    slurmBatchFile = f"{args.output}/slurm.sh"
    slurmSampleList = f"{args.output}/slurmSamples"
    gvcfDir = f"{args.output}/gvcf"

    if (not os.path.isfile(args.sample)):
        parser.print_usage()
        print(f"Error: Sample file {args.sample} does not exist!")
        sys.exit()

    if (not os.path.isdir(args.reference)):
        parser.print_usage()
        print(f"Error: Reference {args.reference} does not exist or is not a directory!")
        sys.exit()

    if (not (os.path.isfile(f"{args.reference}/genome.fasta")  and  os.path.isfile(f"{args.reference}/genome.dict")  and  os.path.isfile(f"{args.reference}/genome.fasta.fai") and  os.path.isfile(f"{args.reference}/genome.fasta.bwt") )):
        parser.print_usage()
        print(f"Error: Reference directory {args.reference} misses some files. Please run script /home/vitisgen/tools/gatk_genome.sh  to create this directory!")
        sys.exit()


    if (not os.path.exists(args.output)):
        os.mkdir(args.output)

    if (not os.path.exists(gvcfDir)):
        os.mkdir(gvcfDir)

    logging.basicConfig(filename=args.output + '/run.log',level=logging.DEBUG)

    logging.debug(f"Command: {commandLine}")


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
            if (not re.search("\w", line)):
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
                    sampleName = sampleName + "__" + plateWell

            sampleList.append(sampleName)

            if ((not os.path.isfile(fieldArray[2])) and  ("1" not in args.skip)):
                print(f"Error: Sample fastq file {fieldArray[2]} does not exist!")
                sys.exit()
            if ((not os.path.isfile(fieldArray[3])) and ("1" not in args.skip)):
                print(f"Error: Sample fastq file {fieldArray[3]} does not exist!")
                sys.exit()
                
            sampleToFileList.append((sampleName, fieldArray[2], fieldArray[3]))

    # run haplotype caller
    if ("1" not in args.skip):
        logging.info("Step 1: run gatk haplotypecaller")
        run_haplotypecaller()

    # get sequence tag list across population
    if ("2" not in args.skip):
        logging.info("Step 2: run genotyper")
        genotyper()


def run_haplotypecaller():
    print("Run run_haplotypecaller ")
    if (args.slurmcluster == ""):
        print("on single node ... ")
        sFh = open (slurmSampleList, "w")
        for ss in sampleToFileList:
            sFh.write(f"{ss[0]}\t{ss[1]}\t{ss[2]}\n")
        sFh.close()
#        print ("Checking dependencies: ")
#        if (not checkApp("gatk")):
#            sys.exit()
#
#        if (not checkApp("bwa")):
#            sys.exit()
#
#        if (not checkApp("samtools")):
#            sys.exit()

        pool = multiprocessing.Pool(processes= args.job)
        pool.starmap(haplotypecaller, sampleToFileList)
        pool.close()


    else:
        print(f"on slurm cluster {args.slurmcluster} ... ")

        ## check finished files for restart to be restartable
        processedList = glob.glob(f"{gvcfDir}/*.done")
        processedList =  list(map(lambda each:each.replace(f"{gvcfDir}/", ""), processedList))
        processedList =  set(map(lambda each:each.replace(".done", ""), processedList))
        print ("Finished samples (will be skipped):")
        print (processedList)

        hostName =os.uname()[1]
        curr_wd = os.getcwd()

        slurmSampleListabs = os.path.abspath(slurmSampleList)
        gvcfDirabs = os.path.abspath(gvcfDir)
        refDirabs = os.path.abspath(args.reference)

        jobCounts =0 
        sFh = open (slurmSampleList, "w")
        for ss in sampleToFileList:
            if (not ss[0] in processedList):
                jobCounts += 1
                sFh.write(f"{ss[0]}\t{ss[1]}\t{ss[2]}\n")
        sFh.close()

        if (jobCounts > 0):
            jobCounts = jobCounts -1

            sFh = open (slurmBatchFile, "w")
            sFh.write(f'''#!/bin/bash
#
#SBATCH --job-name=slurmgatk
#SBATCH --cluster={args.slurmcluster}
#SBATCH --chdir=/workdir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={args.thread}
#SBATCH --time=10-0
#SBATCH --mem-per-cpu=1G
#SBATCH --array=0-{jobCounts}:{args.slurmBatchSize}

srun {slurmScript} {hostName} {curr_wd} {slurmSampleListabs} {refDirabs} {gvcfDirabs} {args.trim} {args.slurmBatchSize} {args.thread} $SLURM_ARRAY_TASK_ID

''')
            sFh.close()
            #os.system(f"sbatch {slurmBatchFile}")
            sys.exit(f"Run slurm command manually: sbatch {args.output}/slurm.sh")         
        else:
            logging.info("cutadapt slurms jobs all finished.")


def haplotypecaller(sampleName, file1, file2):
    workdir=f"{gvcfDir}/{sampleName}"
    pathlib.Path(workdir).mkdir(parents=True, exist_ok=True)

    cmd = f"{bbmerge} t={args.thread} in1={file1} in2={file2} outm={workdir}/contig.fastq"
    bmerge_returned_value = os.system(cmd)  # returns the exit code in unix

    if (args.trim > 0):
        wFh = open (f"{workdir}/tmp.fastq", "w")
        with open(f"{workdir}/contig.fastq", 'r') as fh:
            for line in fh:
                if (line.startswith("@") or line.startswith("+")):
                    wFh.write(line)
                    continue
                line = line.rstrip()
                if ((args.trim *2) < len(line)):
                    line = line[args.trim : len(line)-args.trim]
                wFh.write(line + "\n")
        fh.close()
        wFh.close()
        os.system(f"mv {workdir}/contig.fastq {workdir}/contig.ori.fastq")
        os.system(f"mv {workdir}/tmp.fastq {workdir}/contig.fastq")


    cmd = f"{bwa} mem -v 1 -t {args.thread} {args.reference}/genome.fasta {workdir}/contig.fastq -R '@RG\\tID:{sampleName}\\tSM:{sampleName}' > {workdir}/pass1.sam"
    bwa_returned_value = os.system(cmd)  # returns the exit code in unix

    #filter sam file
    wFh = open (f"{workdir}/pass2.sam", "w")
    with open(f"{workdir}/pass1.sam", 'r') as fh:
        for line in fh:
            if (re.match("@", line)):
                wFh.write(line)
                continue
            fieldArray = line.split(sep="\t")
            flag = fieldArray[1]

            #filter by flag
            if (not ((flag == "0" ) or (flag == "16" ))):
                continue

            # filter by cigar soft clipping
            cigar =  fieldArray[5]
            match = re.findall("(\d+)S", cigar)
            if (len(match) == 0):
                pass
            else:
                #maxClipLen = max(list(map(int, match)))
                sumClipLen = sum(list(map(int, match)))
                if (sumClipLen >5):
                    continue


            #filter by indel
            indelCount = len(re.findall("[ID]", cigar))
            if (indelCount>2):
                continue

            # filter by alignment score normalized by read length 
            match = re.search("AS:i:(\d+)", fieldArray[13])
            if (match):
                algnCore = int(match[1])
                seqlen = len(fieldArray[9])
                normalizedScore = algnCore/seqlen
                if (normalizedScore<0.7):
                    #print(f"{algnCore} {fieldArray[11]} {cigar} {normalizedScore}")
                    continue 
                #elif (normalizedScore<0.8):
                #    print(f"{algnCore} {seqlen} {fieldArray[11]} {cigar} {normalizedScore}")
            else:
                continue  

            #filter by mapq
            mapq = fieldArray[4]
            if (int(mapq) <40):
                continue


            wFh.write(line)
    fh.close()
    wFh.close()

    cmd = f"{samtools} sort -T {workdir} -m 1G -O bam -o {workdir}/pass3.bam {workdir}/pass2.sam \n"
    cmd += f"{picard} CleanSam INPUT={workdir}/pass3.bam OUTPUT={workdir}/pass4.bam \n"
    cmd += f"mv {workdir}/pass4.bam {gvcfDir}/{sampleName}.bam \n"
    cmd += f"{picard} BuildBamIndex INPUT={gvcfDir}/{sampleName}.bam \n"

    returned_value = os.system(cmd)
    if (returned_value==0):
        logging.info("clean bam done: " + sampleName)
    else:
        logging.warning("clean bam error: " + sampleName)

    cmd = f"{gatk} --java-options \"-Djava.io.tmpdir={gvcfDir}\"  HaplotypeCaller -R {args.reference}/genome.fasta -I {gvcfDir}/{sampleName}.bam -ERC GVCF --native-pair-hmm-threads {args.thread} -O {gvcfDir}/{sampleName}.g.vcf"
    logging.info("gatk cmd: " +cmd)
    returned_value = os.system(cmd)

    os.system(f"rm -fr {workdir}")

def genotyper():

    print ("Checking dependencies: ")
    #if (not checkApp("gatk")):
    #    sys.exit()


    sampleListStr = ""
    for sampleName in sampleList:
        sampleListStr += f" --variant {gvcfDir}/{sampleName}.g.vcf"
    cmd = f"{gatk} --java-options \"-Djava.io.tmpdir={gvcfDir} -Xmx8g\"  CombineGVCFs -R {args.reference}/genome.fasta {sampleListStr} -O {gvcfDir}/all.g.vcf -L chr9|19496735-19496988"
    returned_value = subprocess.call(cmd, shell=True)
    cmd = f"{gatk} --java-options \"-Djava.io.tmpdir={gvcfDir} -Xmx8g\"  GenotypeGVCFs  -R {args.reference}/genome.fasta -V {gvcfDir}/all.g.vcf -O {gvcfDir}/gatk.vcf -L chr9|19496735-19496988"
    returned_value = subprocess.call(cmd, shell=True)

    os.system ("rm {gvcfDir}/*.bam*")
    os.system ("rm {gvcfDir}/*.g.vcf*")


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


