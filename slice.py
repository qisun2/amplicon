#!/usr/bin/env python3
import sys
import argparse
import os
from os.path import isfile, join
import re

def main():
    if sys.version_info[0] < 3:
        raise Exception("This code requires Python 3.")

    parser = argparse.ArgumentParser(description='Run slicer.')
    parser.add_argument('-i','--input',type=str,required=True,help='Input directory, it must contain at least one file named hap_genotype from the amplicon.py tool.')
    parser.add_argument('-o','--output',type=str,required=True,help='Output directory.')
    parser.add_argument('-f','--familyFile',type=str,required=True,help='Family file. It is a text file, with one individual per line. The individual name must be in format plateName_well.')
    #parser.add_argument('-m','--familyName',type=str,required=False,help='Family name. A string with no space.')
    parser.add_argument('-a','--alleleFreq',type=float,required=True,help='Minimum allele frequency.')

    #global variables
    global args
    global sampleNameDupCheck
    global plateWellList
    global plateWellToSample
    global sampleToPlateWell
    global polidy 
    global matrixCountProcessed
    global hapGenotypeProcessed

    matrixCountProcessed=False
    hapGenotypeProcessed=False
    
    polidy=2

    args=parser.parse_args()

    if (not os.path.isfile(args.familyFile)):
        parser.print_usage()
        print(f"Error: family file {args.familyFile} does not exist!")
        sys.exit()


    #process family file
    plateWellList = []
    plateWellToSample = {}
    sampleToPlateWell = {}

    with open(args.familyFile, 'r') as fhs:
        for line in fhs:
            if (re.search("\w", line)):
                line = line.strip()
                df = line.split("\t")
                line =df[0]
                if line not in plateWellList:
                    plateWellList.append(line)
                continue
        fhs.close()
  


    ### create output directories
    try:
        os.mkdir(args.output)
    except OSError:
        print ("Creation of the output directory %s already exists. Please delete the directory." % args.output)
        sys.exit()
    else:
        print ("Output directory created: %s " % args.output)


    #process countMatrix file
    readCountMatrixFile = f"{args.input}/markerToSampleReadCountMatrix"
    if os.path.isfile(readCountMatrixFile):
        fhs =  open(readCountMatrixFile, 'r')
        ### read header 
        MarkerLine = fhs.readline()

        filePath =  f"{args.output}/readCountMatrixFile"
        fh = open (filePath,'w')
        fh.write(MarkerLine)
        
        for line in fhs:
            if (not re.search("\w", line)):
                continue
            fieldArray = line.split("\t", 1)
            sampleName = fieldArray[0]
            if ("__" not in sampleName):
                continue
            plateWell = sampleName.split("__")[1]

            if plateWell not in plateWellList:
                continue
            
            plateWellToSample[plateWell] = sampleName
            sampleToPlateWell[sampleName] = plateWell
            fh.write(line)  
        fhs.close()
        fh.close()

        t1= len(plateWellList)
        print (f"Number of individuals in List: {t1}")
        t2 = len(plateWellToSample)
        print (f"Number of individuals found in data: {t2}")

        if (t2 < t1):
            print (f"The following individuals are not found in the amplicon.py output directory {args.input}. Please correct them and try again:")
            for t in plateWellList:
                if t not in plateWellToSample:
                    print(t)
            sys.exit()
        matrixCountProcessed=True
        
        
        
    #process hap_genotype file    

    hap_genotype = f"{args.input}/hap_genotype"
    if (not os.path.isfile(hap_genotype)):
        parser.print_usage()
        print(f"Error: Read matrix file {hap_genotype} does not exist!")
        sys.exit()

    fhs =  open(hap_genotype, 'r')


    ### read header 
    sampleLine = fhs.readline()
    sampleLine = sampleLine.rstrip()
    sampleList = sampleLine.split(sep="\t")
    tt = sampleList.pop(0)
    tt = sampleList.pop(0)

    sampleNameToIndex = {}
    plateWellToSample={}
    sampleToPlateWell={}
    indexList = []
    sIndex =0

    for sampleName in sampleList:
    
        if ("__" not in sampleName):
            continue
        plateWell = sampleName.split("__")[1]
            
            
        if plateWell in plateWellList:
            sampleNameToIndex[sampleName] = sIndex
            plateWellToSample[plateWell] = sampleName
            sampleToPlateWell[sampleName] = plateWell
        
        sIndex+=1
    
    #first make sure the sample list provide in familyFile are all present in the hap_genotype file
    t1= len(plateWellList)
    print (f"Number of individuals in List: {t1}")
    t2 = len(plateWellToSample)
    print (f"Number of individuals found in data: {t2}")

    if (t2 < t1):
        print (f"The following individuals are not found in the amplicon.py output directory {args.input}. Please correct them and try again:")
        for t in plateWellList:
            if t not in plateWellToSample:
                print(t)
        sys.exit()
        
    filePath =  f"{args.output}/hap_genotype"
    fh = open (filePath,'w')
    fh.write ("Locus\tHaplotypes")
    
    for plateWell in plateWellList:
        sampleName = plateWellToSample[plateWell]
        sIndex = sampleNameToIndex[sampleName]
        indexList.append(sIndex)
        fh.write(f"\t{sampleName}")
    fh.write(f"\n")

    for line in fhs:
        if (not re.search("\w", line)):
            continue
        line = line.rstrip("\n")
        fieldArray = line.split(sep="\t")
        markerName = fieldArray.pop(0)
        haps = fieldArray.pop(0)
        totalGametes = 0
        alleleCounts = {}
        keptAlleles = {}

        ## calculate allele frequency for each allele, and determine used allele for family based on frequency
        for i in indexList:
            GT = fieldArray[i].split(":")                
            if (GT[0] == "./."):
                pass
            else:
                totalGametes += 2
                alleles = GT[0].split("/")
                for a in alleles:
                    if a in alleleCounts:
                        alleleCounts[a] +=1
                    else:
                        alleleCounts[a] =1
        if (totalGametes==0):
            print(f"{markerName}: no data")

        sortedAlleles = sorted(alleleCounts.items(), key=lambda kv: kv[1], reverse=True)
        newHapStr = "";
        newTotalGametes = 0
        for alleleCount in sortedAlleles:
            aId = alleleCount[0]
            frq = int(alleleCount[1])/totalGametes
            if (frq > args.alleleFreq):
                keptAlleles[aId] = 1
                newTotalGametes += int(alleleCount[1])
        
        if (newTotalGametes == 0):
            print(f"{markerName}: no used alleles")

        # re-calculate allele frq based on used allelels
        for alleleCount in sortedAlleles:
            aId = alleleCount[0]
            if aId in keptAlleles:
                b = "%4.2f" %(int(alleleCount[1])/newTotalGametes)
                newHapStr += f"{aId}({b});"  
        fh.write(f"{markerName}\t{newHapStr}")

        # rewrite genotype based on used alleles
        GTstr = ""
        for i in indexList:
            GT = fieldArray[i].split(":")                
            if (GT[0] == "./."):
                GTstr = "./.:0"
            else:
                alleles = GT[0].split("/")
                if (alleles[0] == alleles[1]):
                    if alleles[0] in keptAlleles:
                        GTstr=fieldArray[i]
                    else:
                        GTstr = "./.:0"
                else:
                    readCounts = GT[1].split(",")
                    kept = []
                    keptCounts = []
                    for a in [0,1]:
                        if alleles[a] in keptAlleles:
                            kept.append(alleles[a])
                            keptCounts.append(readCounts[a])
                    if (len(kept) == 0):
                        GTstr = "./.:0"
                    elif (len(kept) == 1):
                        GTstr=f"{kept[0]}/{kept[0]}:{keptCounts[0]}"
                    else:
                        GTstr=f"{kept[0]}/{kept[1]}:{keptCounts[0]},{keptCounts[1]}"
            fh.write(f"\t{GTstr}")
        fh.write("\n")
    fhs.close()
    fh.close()

if __name__=="__main__":
    main()


