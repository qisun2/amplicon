#!/usr/bin/env python3

import os
import sys
import pathlib
import shutil
import subprocess
import re

bwa = "bwa"
bbmerge = "/programs/bbmap-38.45/bbmerge.sh"
gatk = "/programs/gatk4/gatk"
samtools = "samtools"
picard = "java -jar /programs/picard-tools-2.19.2/picard.jar"
baseWorkdir = "/workdir"
RELEASE = "/programs/sentieon-genomics-201808.05"
os.environ["SENTIEON_LICENSE"] = "cbsulogin2.tc.cornell.edu:8990"

def main():
    global hostName
    global datadir
    global markerList
    global destResultDir
    global workdir
    global sampleFile
    global threads
    global trim
    global slurmBatchSize
    global refDir

    hostName=sys.argv[1]
    datadir=sys.argv[2]
    inputSampleFile = sys.argv[3]
    inputRefDir = sys.argv[4]
    destResultDir = sys.argv[5]
    trim = int(sys.argv[6])
    slurmBatchSize = int(sys.argv[7])
    threads = int(sys.argv[8])
    sampleID = int(sys.argv[9])

    workdir= baseWorkdir + "/ampseq" + str(sampleID)
    sampleFile =  workdir + "/sampleFile"
    refDir = workdir + "/ref"

    os.system(f"rm -rf {workdir}")
    pathlib.Path(workdir).mkdir(parents=True, exist_ok=True)
    os.chdir(workdir)

    # scp the sample list file
    while (True):
        return_value = os.system(f"scp -q -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null {hostName}:{inputSampleFile} {sampleFile}")
        if (return_value==0):
            break
        time.sleep(30)
 
    # scp the reference file
    while (True):
        return_value = os.system(f"scp -r -q -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null {hostName}:{inputRefDir} {refDir}")
        if (return_value==0):
            break
        time.sleep(30)


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

        haplotypecaller(sampleName, file1, file2)

    os.chdir("../")
    os.system(f"rm -rf {workdir}")


def haplotypecaller(sampleName, file1, file2):

    while (True):
        cmd = f"scp -q -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null {hostName}:{datadir}/{file1} {workdir}/ "
        scp_returned_value = os.system(cmd)
        if (scp_returned_value==0):
            print(f"scp return {scp_returned_value} for {sampleName}", flush=True)
            break
        print(f"scp return {scp_returned_value} for {sampleName}; wait 60 seconds and try again", flush=True)
        time.sleep(60)

    while (True):
        cmd = f"scp -q -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null {hostName}:{datadir}/{file2} {workdir}/ "
        scp_returned_value = os.system(cmd)
        if (scp_returned_value==0):
            print(f"scp return {scp_returned_value} for {sampleName}", flush=True)
            break
        print(f"scp return {scp_returned_value} for {sampleName}; wait 60 seconds and try again", flush=True)
        time.sleep(60)


    cmd = f"{bbmerge} t={threads} in1={file1} in2={file2} outm={workdir}/contig.fastq"
    bmerge_returned_value = os.system(cmd)  # returns the exit code in unix


    if (trim > 0):
        wFh = open (f"{workdir}/tmp.fastq", "w")
        with open(f"{workdir}/contig.fastq", 'r') as fh:
            for line in fh:
                if (line.startswith("@") or line.startswith("+")):
                    wFh.write(line)
                    continue
                line = line.rstrip()
                if ((trim *2) < len(line)):
                    line = line[trim : len(line)-trim]
                wFh.write(line + "\n")
        fh.close()
        wFh.close()
        os.system(f"mv {workdir}/contig.fastq {workdir}/contig.ori.fastq")
        os.system(f"mv {workdir}/tmp.fastq {workdir}/contig.fastq")


    cmd = f"{bwa} mem -v 1 -t {threads} {refDir}/genome.fasta {workdir}/contig.fastq -R '@RG\\tID:{sampleName}\\tSM:{sampleName}' > {workdir}/pass1.sam"
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
    cmd += f"mv {workdir}/pass4.bam {workdir}/{sampleName}.bam \n"
    cmd += f"{picard} BuildBamIndex INPUT={workdir}/{sampleName}.bam \n"

    clean_returned_value = os.system(cmd)

    #cmd = f"{gatk} --java-options \"-Djava.io.tmpdir={workdir}\"  HaplotypeCaller -R {refDir}/genome.fasta -I {workdir}/{sampleName}.bam -ERC GVCF --native-pair-hmm-threads {threads} -O {workdir}/{sampleName}.g.vcf"
    cmd = f"{RELEASE}/bin/sentieon driver -t {threads}  -r {refDir}/genome.fasta -i {workdir}/{sampleName}.bam --algo Haplotyper --min_map_qual 50 --emit_mode gvcf {workdir}/{sampleName}.g.vcf";
    gatk_returned_value = os.system(cmd)

    if ((scp_returned_value==0) and (bwa_returned_value == 0) and (clean_returned_value == 0) and (gatk_returned_value == 0)):
        cmd = f"scp -q -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null {workdir}/{sampleName}.g.vcf* {hostName}:{destResultDir} \n"
        cmd = cmd + f"touch ./{sampleName}.done  \n"
        cmd = cmd + f"scp -q -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null ./{sampleName}.done {hostName}:{destResultDir} \n"
        os.system(cmd)

########### end of  PCR error correction code #################################

if __name__=="__main__":
    main()


