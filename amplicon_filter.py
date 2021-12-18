#!/usr/bin/env python3


from os.path import isfile, join, exists
import argparse
import os
import sys
import subprocess
import re
import swalign

from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import SeqIO
from collections import defaultdict



blastCmd = "blastn"
makeblastdbCmd = "makeblastdb"
tmpDir = "./"
refFilterLength = 25
def main():

    global args
    global markerList
    global refSeq
    global alleleSeqs
    global acceptedAlleleIds

    commandLine = " ".join(sys.argv[1:])

    parser = argparse.ArgumentParser(description='Run GATK Haplotype Caller on ampseq data.')

    # Required arguments
    parser.add_argument('-a','--alleleFasta',type=str,required=True,help='Allele fasta file. It should be the output file HaplotypeAllele.fasta from amplicon.py.')
    parser.add_argument('-r','--refSeq',type=str,required=False,default="",help='Reference sequence for each marker in fasta format. If the value is not set, the first allele (most common allele) is used as refseq')
    parser.add_argument('-o','--outputFasta',type=str,required=True,help='Output fasta file. It is required')
    parser.add_argument('-f','--filter',type=str,required=False,default="1",help='Filter alleles by reference sequence. 1. No filter; 2 Filter by smith-waterman alignment to reference; 3. Filter by BLAST. default 1')
    parser.add_argument('-x','--alnpct',type=float,required=False,default=0.7,help='percent identity of alignment')
    parser.add_argument('-e','--evalue',type=float,required=False,default=1e-2,help='Blast to reference sequence evalue cutoff')
    parser.add_argument('-m','--removed',type=str,required=False,default=f"{tmpDir}/amplicon.removed.fasta",help='removed sequences')


    if sys.version_info[0] < 3:
        raise Exception("This code requires Python 3.")

    args=parser.parse_args()
 

    # check reference sequence file exist if parameter set
    if (args.filter not in ["1", "2", "3"]):
        print (f"Error: The --filter option must be a value of 1, 2 or 3")
        parser.print_usage()
        sys.exit()  

    #validate input alleleFasta file
    alleleFasta = args.alleleFasta
    if (not os.path.isfile(alleleFasta)):
        print(f"Error: Reference sequence fasta file {refseqFile} does not contain reference sequence for marker {m}!")
        sys.exit()
    
    #get marker List, and retrieve first allele for each marker as reference
    markerList = []
    refSeqs = {}
    alleleSeqs = {}
    for record in SeqIO.parse(alleleFasta, "fasta"):
        matches = re.match("(\\S+)#(\\S+)", record.id)
        if (matches):
            markerId = matches[1]
            alleleId = matches[2]
            if (markerId not in refSeqs):
                refSeqs[markerId] = record
                markerList.append(markerId)
                alleleSeqs[markerId] = []
            alleleSeqs[markerId].append(record)
        else:
            print(f"Error: the fasta file is not an output from amplicon.py. The sequence id {record.id} does not match the expected format!")
            sys.exit()


    #if reference file  exists, set reference sequence from refseq file
    refseqFile = args.refSeq
    t = set()
    if (refseqFile!="") and (os.path.isfile(refseqFile)):
        refSeqs = {}
        for record in SeqIO.parse(refseqFile, "fasta"):
            refSeqs[record.id] = record
        for m in markerList:
            if m not in refSeqs:
                print(f"Error: Reference sequence fasta file {refseqFile} does not contain reference sequence for marker {m}!")
                sys.exit()


    acceptedAlleleIds = set()
    if (args.filter=="1"):
        print (f"Filtering method is set to no filtering. For filtering, -f should have a value 2(Smith Waterman) or 3 (BLAST)")
        sys.exit()
    elif (args.filter=="2" or args.filter=="3"):
        for markerId in markerList:
            refSeqStr = str(refSeqs[markerId])
            seqLen = len(refSeqStr)
            refSeq_5p25nt = seqStr[0:refFilterLength]
            refSeq_3p25nt = seqStr[seqLen-refFilterLength:seqLen]
            if args.filter=="2":
                filterBySwalign(markerId, refSeq_5p25nt, refSeq_3p25nt)
            else:
                filterByBlast(markerId, refSeq_5p25nt, refSeq_3p25nt)
    else:
        print (f"{args.filter} is not valid for filtering.  -f should have a value 2(Smith Waterman) or 3 (BLAST)")
        parser.print_usage()
        sys.exit()

    # write to output
    wh = open(outputFasta, "w")
    rwh = open(args.removed, "w")
    for markerId in markerList:
        for record in alleleSeqs[markerId]:
            if record.id in acceptedAlleleIds:
                SeqIO.write(record, wh, "fasta")
            else:
                SeqIO.write(record, rwh, "fasta")
    wh.close()
    rwh.close()



def filterBySwalign(markerId, refSeq_5p25nt, refSeq_3p25nt):
    matchedNT = args.alnpct*refFilterLength
    alleleRecords = alleleSeqs[markerId]
    for record in alleleRecords:
        matches = re.match("(\\S+)#(\\S+)", record.id)
        seqStr = str(record.seq)
        seqLen = len(seqStr)
        query_5p25nt = seqStr[0:refFilterLength]
        query_3p25nt = seqStr[seqLen-refFilterLength:seqLen]
        alignment1 = sw.align(refseq_5p25nt, query_5p25nt)
        alignment2 = sw.align(refseq_3p25nt, query_3p25nt)
        if (alignment1.matches>matchedNT and alignment2.matches>matchedNT):
            acceptedAlleleIds.add(record.id)
       

def filterByBlast(markerId, refSeq_5p25nt, refSeq_3p25nt):
        ref5Query = f"{tmpDir}/tmp.ref5.fasta"
        ref3Query = f"{tmpDir}/tmp.ref3.fasta"
        blastOut5 = f"{tmpDir}/tmp.ref5.blast"
        blastOut3 = f"{tmpDir}/tmp.ref3.blast"
        dbFasta = f"{tmpDir}/tmp.db.fasta"

        w1 = open(ref5Query, "w")
        w1.write(f">ref5\n{refseq_5p25nt}\n")
        w1.close()

        w1 = open(ref3Query, "w")
        w1.write(f">ref3\n{refseq_3p25nt}\n")
        w1.close()

        w1 = open(dbFasta, "w")
        SeqIO.write(alleleSeqs[markerId], w1, "fasta")
        w1.close()

        cmd = f"{makeblastdbCmd} -in {dbFasta} -dbtype nucl"
        subprocess.call(cmd, shell=True)
        cmd = f"{blastCmd} -query {ref5Query} -db {dbFasta} -task \"blastn-short\" -outfmt 6 -out {blastOut5} -max_target_seqs 2000"
        subprocess.call(cmd, shell=True)
        cmd = f"{blastCmd} -query {ref3Query} -db {dbFasta} -task \"blastn-short\" -outfmt 6 -out {blastOut3} -max_target_seqs 2000"
        subprocess.call(cmd, shell=True)

        hits5 = []
        hits3 = []
        with open (blastOut5, "r") as BLASTIN:
            for line in BLASTIN:
                dataFields = line.split("\t")
                if (float(dataFields[2])/100 >= alignpct ):
                    hits5.append(dataFields[1])
            BLASTIN.close()
        with open (blastOut3, "r") as BLASTIN:
            for line in BLASTIN:
                dataFields = line.split("\t")
                if (float(dataFields[2])/100 >= alignpct ):
                    hits3.append(dataFields[1])
            BLASTIN.close()

        hits = [value for value in hits5 if value in hits3]
        for h in hits:
            acceptedAlleleIds.add(h)



if __name__=="__main__":
    main()


