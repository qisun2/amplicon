#!/usr/bin/env python3

import os
import sys
import pathlib
import shutil
import subprocess
import re


bbmergeCMD = "/programs/bbmap-38.45/bbmerge.sh"
cutadaptCMD = "cutadapt"
baseWorkdir = "/workdir"


def main():
    global hostName
    global datadir
    global markerList
    global destResultDir
    global minHaplotypeLength
    global maxHaplotypePerSample
    global maxAlleleReadCountRatio
    global workdir
    global sampleFile
    global primerFile
    global threads
    global slurmBatchSize

    
    hostName=sys.argv[1]
    datadir=sys.argv[2]
    inputSampleFile = sys.argv[3]
    inputPrimerFile = sys.argv[4]
    destResultDir = sys.argv[5]
    minHaplotypeLength = sys.argv[6]
    maxHaplotypePerSample = int(sys.argv[7])
    maxAlleleReadCountRatio = int(sys.argv[8])
    slurmBatchSize = int(sys.argv[9])
    threads = int(sys.argv[10])
    sampleID = int(sys.argv[11])

    workdir= baseWorkdir + "/ampseq" + str(sampleID)
    sampleFile =  workdir + "/sampleFile"
    primerFile = workdir + "/primer.fa"


    markerList = []

    os.system(f"rm -rf {workdir}")
    pathlib.Path(workdir).mkdir(parents=True, exist_ok=True)
    os.chdir(workdir)

    # scp the sample list file
    os.system(f"scp -q -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null {hostName}:{inputSampleFile} {sampleFile}")
    os.system(f"scp -q -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null {hostName}:{inputPrimerFile} {primerFile}")

    with open(primerFile, 'r') as fhs:
        for line in fhs:            
            if (line.startswith( '>' )):
                line = line.rstrip()
                markerList.append(line.replace(">", ""))

    sampleList = []

    with open(sampleFile, 'r') as fhs:
        for line in fhs:
            if (not re.match("\w", line)):
                continue
            line = line.rstrip()
            fieldArray = line.split(sep="\t")
            sampleName = re.sub("\s", "", fieldArray[0])
            file1 = re.sub("\s", "", fieldArray[1])
            file2 = re.sub("\s", "", fieldArray[2])
            sampleList.append((sampleName, file1, file2))

    upperIndex  = sampleID+slurmBatchSize
    if ( upperIndex> len(sampleList)):
        upperIndex = len(sampleList)
    for i in range (sampleID, upperIndex):
        sss = sampleList[i]
        sampleName = sss[0]
        file1 = sss[1]
        file2 = sss[2]

        splitByCutadapt(sampleName, file1, file2)

    os.chdir("../")
    os.system(f"rm -rf {workdir}")



def splitByCutadapt(sampleName, file1, file2):
    cmd = f"scp -q -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null {hostName}:{datadir}/{file1} {workdir}/ \n"
    cmd += f"scp -q -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null {hostName}:{datadir}/{file2} {workdir}/ \n"
    scp_returned_value = os.system(cmd)
    print(f"scp return {scp_returned_value} for {sampleName}", flush=True)

    try:
        sampleDir = f"./{sampleName}"
        if (not os.path.exists(sampleDir)):
            os.mkdir(sampleDir)

        #contig the paired end reads
        cmd = f"{bbmergeCMD} t={threads} in1={file1} in2={file2} outm={sampleDir}/contig.fastq"
        #returned_value_bbmerge = subprocess.call(cmd, shell=True)
        returned_value_bbmerge = os.system(cmd)
        print(f"bbmerge return {returned_value_bbmerge} for {sampleName}", flush=True)

        #run cutadapt to demultiplexing by primers
        cmd = f"{cutadaptCMD} --quiet -e 0.1 --minimum-length={minHaplotypeLength} --trimmed-only -g file:{primerFile} -o {sampleDir}/{{name}}.fastq {sampleDir}/contig.fastq "
        #returned_value_cutadapt = subprocess.call(cmd, shell=True)
        returned_value_cutadapt = os.system(cmd)
        print(f"cutadapt return {returned_value_cutadapt} for {sampleName}", flush=True)

        #collapse identical reads
        tagBySampleFile = f"./{sampleName}.tbs"
        tbsfh = open(tagBySampleFile, "w")

        readCountMatrix = {}
        
        rcFh = open (f"./{sampleName}.readcount", "w")
        rcFh.write(f"{sampleName}")
        for marker in markerList:
            readCount =0
            markerFile = f"{sampleDir}/{marker}.fastq"
            collapsedMarkerFile = f"{sampleDir}/{marker}.collapsed"
            if (os.path.exists(markerFile)):
                cmd=f"awk 'NR%4==2' {markerFile}|LC_ALL=C sort |uniq -c |LC_ALL=C  sort -k1,1rn > {collapsedMarkerFile}"
                #returned_value = subprocess.call(cmd, shell=True)
                awk_returned_value = os.system(cmd)
                #print(f"collapse return {awk_returned_value} for {sampleName}")
                tagCount = 0
                topAlleleReadCount = 0
                
                with open(collapsedMarkerFile, 'r') as cfh:
                    for line in cfh:
                        [copyNumber, seqStr] = line.split()
                        copyNumber = int(copyNumber)
                        readCount += copyNumber
                        if (tagCount==0):
                            topAlleleReadCount = copyNumber
                        if (tagCount < maxHaplotypePerSample):
                            if ((topAlleleReadCount/copyNumber) < maxAlleleReadCountRatio):
                                tagCount+=1
                                tbsfh.write(f"{sampleName}\t{marker}\t{seqStr}\t{copyNumber}\n")
                    cfh.close()
            rcFh.write(f"\t{readCount}")
        tbsfh.close()
        rcFh.write("\n")
        rcFh.close()


        if ((scp_returned_value==0) and (returned_value_bbmerge == 0) and (returned_value_cutadapt == 0)):
            cmd = f"scp -q -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null {tagBySampleFile} {hostName}:{destResultDir} \n"
            cmd = cmd + f"scp -q -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null ./{sampleName}.readcount {hostName}:{destResultDir} \n"
            cmd = cmd + f"touch ./{sampleName}.done  \n"
            cmd = cmd + f"scp -q -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null ./{sampleName}.done {hostName}:{destResultDir} \n"
            os.system(cmd)
        return 1 
    except Exception as e:
        print(f'Caught exception in job {sampleName} {marker}')
        traceback.print_exc()
        print()
        raise e



########### end of  PCR error correction code #################################

if __name__=="__main__":
    main()


