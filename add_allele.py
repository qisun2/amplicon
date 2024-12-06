#!/usr/bin/env python3


from Bio import SeqIO
import re
import argparse

parser = argparse.ArgumentParser(description='add new alleles')

# Required arguments
parser.add_argument('-o','--oldfile',type=str,required=True,help='old haplotypeallele.fasta file')
 
parser.add_argument('-n','--newfile',type=str,required=True,help='new haplotypeallele.fasta file')
    
parser.add_argument('-m','--mergedfile',type=str,required=True,help='merged haplotypeallele.fasta file')

args=parser.parse_args()

#import original file
oriFile = args.oldfile
newFile = args.newfile  #new output file from amplicon.py
newHapFile = args.mergedfile #a new merged haplotypeallele.fasta file
statFile = "merging.stat" #stats of new alleles

marker2allele = {}
marker2MaxAlleleID = {} 
marker2newAlleleCount = {} 
marker2alleleForNewFile = {}

allMarkers = []

with open(oriFile) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        marker, alleleId = record.id.split("#")
        alleleId = int(alleleId)
        #print(marker, alleleId)
        if marker not in marker2allele:
            allMarkers.append(marker)
            marker2allele[marker] = {}
            marker2alleleForNewFile[marker] = ""
            marker2newAlleleCount[marker] =0
            marker2MaxAlleleID[marker]=alleleId
        marker2allele[marker][record.seq] = alleleId
        marker2alleleForNewFile[marker] += f">{marker}#{alleleId}\n{record.seq}\n"
        if marker2MaxAlleleID[marker] < alleleId:
            marker2MaxAlleleID[marker] = alleleId
        

with open(newFile) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        marker, alleleId = record.id.split("#")
        #print(marker, alleleId)
        if marker not in marker2allele:
            allMarkers.append(marker)
            marker2allele[marker]={}
            marker2MaxAlleleID[marker] = 0
            marker2newAlleleCount[marker] = 0
            marker2alleleForNewFile[marker] = ""
        if record.seq not in marker2allele[marker]:
            marker2newAlleleCount[marker] += 1
            marker2MaxAlleleID[marker] +=1
            marker2alleleForNewFile[marker] += f">{marker}#{marker2MaxAlleleID[marker]}\n{record.seq}\n"
SH = open (statFile, "wt")
WH = open (newHapFile, "wt")       
for marker in allMarkers:
    oricount = len(marker2allele[marker])
    SH.write(f"{marker}\t{oricount}\t{marker2newAlleleCount[marker]}\n")
    WH.write(marker2alleleForNewFile[marker])
SH.close()
WH.close()
