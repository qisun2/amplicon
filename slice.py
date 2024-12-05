#!/usr/bin/env python3
import sys
import argparse
import os
import gzip
from os.path import isfile, join
import re

def main():
    if sys.version_info[0] < 3:
        raise Exception("This code requires Python 3.")

    parser = argparse.ArgumentParser(description='Run slicer.')
    parser.add_argument('-i','--input',type=str,required=True,help='Input directory, it must contain a file named hap_genotype or hap_genotype.gz from the amplicon.py tool.')
    parser.add_argument('-o','--output',type=str,required=True,help='Output directory.')
    parser.add_argument('-s','--sampleFile',type=str,required=False,help='Sample file. It is a text file, with one individual per line. The individual name must match the hap_genotype file.')
    parser.add_argument('-f','--familyFile',type=str,required=False,help='Family file. This option is for the Vitisgen project, which use the sample name convention sample__plateName_well. It is a text file, with one individual per line. The individual name must be in format plateName_well.')
    parser.add_argument('-p','--plateList',type=str,required=False,help='Plate list file. This option is for the Vitisgen project, which use the sample name convention sample__plateName_well. It is a text file, with one plate per line.')
    parser.add_argument('-a','--alleleFreq',type=float,required=True,help='Minimum allele frequency.')

    #global variables
    global args
    global plateWellList
    global plateList
    global plateWellToSample
    global mode
    global sampleNameList

    

    args=parser.parse_args()

    if ((args.plateList  == None) and (args.familyFile == None) and (args.sampleFile == None)):
        print(f"Error: the --plateList(-f) or --familyFile(-p) or --sampleFile (-s) parameter must be provided!")
        sys.exit()
 
    if ((args.plateList  != None) and (args.familyFile != None)):
        print(f"Error: the --plateList(-f) and --familyFile(-p) cannot be both provided!")
        sys.exit()
        
    if (args.familyFile != None) and (not os.path.isfile(args.familyFile)):
        parser.print_usage()
        print(f"Error: family file {args.familyFile} does not exist!")
        sys.exit()

    if (args.plateList != None) and (not os.path.isfile(args.plateList)):
        parser.print_usage()
        print(f"Error: plate list file {args.plateList} does not exist!")
        sys.exit()
        
    if (args.sampleFile != None) and (not os.path.isfile(args.sampleFile)) :
        parser.print_usage()
        print(f"Error: sample file {args.sampleFile} does not exist!")
        sys.exit()

    #process family file
    plateWellList = []
    plateWellToSample = {} 
    if (args.familyFile != None):
        mode="sample"
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


            
            
    #process plate file
    plateList = []
    if (args.plateList != None):
        mode="plate"
        with open(args.plateList, 'r') as fhs:
            for line in fhs:
                if (re.search("\w", line)):
                    line = line.strip()
                    df = line.split("\t")
                    line =df[0]
                    if line not in plateList:
                        plateList.append(line)
                    continue
            fhs.close() 

    sampleNameList = []
    if (args.sampleFile != None):
        mode="generic_sample"
        with open(args.sampleFile, 'r') as fhs:
            for line in fhs:
                if (re.search("\w", line)):
                    line = line.strip()
                    df = line.split("\t")
                    line =df[0]
                    if line not in sampleNameList:
                        sampleNameList.append(line)
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

        filePath =  f"{args.output}/markerToSampleReadCountMatrix"
        fh = open (filePath,'w')
        fh.write(MarkerLine)
        foundPlateDict = {}
        for line in fhs:
            if (not re.search("\w", line)):
                continue
            fieldArray = line.split("\t", 1)
            sampleName = fieldArray[0]
            
            if mode == "generic_sample":
                if (sampleName not in sampleNameList):
                    plateWellToSample[sampleName] = sampleName
                    continue
            else:
                if ("__" not in sampleName):
                    continue
                plateWell = sampleName.split("__")[1]
                m= re.match("(.+)_\w{3}$", plateWell)
                if m:
                    plate = m[1]
                else:
                    print(f"Warning: sample platewell {plateWell} is skipped")
                    continue
                
                if (mode=="sample") and (plateWell not in plateWellList):
                    continue
                
                if (mode=="plate") and (plate not in plateList):
                    continue
                
                if (mode=="plate"):
                    if plate in foundPlateDict:
                        foundPlateDict[plate]+=1
                    else:
                        foundPlateDict[plate]=1
            
                plateWellToSample[plateWell] = sampleName            
            fh.write(line)  
        fhs.close()
        fh.close()

        if (mode=="generic_sample"):
            t1= len(sampleNameList)
            print (f"Number of individuals in List: {t1}")
            t2 = len(plateWellToSample)
            print (f"Number of individuals found in markerToSampleReadCountMatrix: {t2}")

                
        if (mode=="sample"):
            t1= len(sampleNameList)
            print (f"Number of individuals in List: {t1}")
            t2 = len(plateWellToSample)
            print (f"Number of individuals found in data: {t2}")

            if (t2 < t1):
                print (f"The following individuals are not found in the amplicon.py output directory {args.input}. Please correct them and try again:")
                for t in sampleNameList:
                    if t not in plateWellToSample:
                        print(t)
                sys.exit()
           
        if (mode=="plate"):
            print ("Found samples per plate in markerToSampleReadCountMatrix:")
            for p in plateList:
                if p not in foundPlateDict:
                    print(f"{p}\t0")
                else:
                    print(f"{p}\t{foundPlateDict[p]}")
      
    #process hap_genotype file    

    hap_genotype = f"{args.input}/hap_genotype"
    hap_genotype_gz = f"{args.input}/hap_genotype.gz"
    
    if (not os.path.isfile(hap_genotype)) and (not os.path.isfile(hap_genotype_gz)):
        parser.print_usage()
        print(f"Error: Read matrix file {hap_genotype} does not exist!")
        sys.exit()
    has_gz = False
    if os.path.isfile(hap_genotype_gz):
        has_gz = True
        
    
    if has_gz:
        fhs =  gzip.open(hap_genotype_gz, 'rt')
    else:
        fhs =  open(hap_genotype, 'rt')


    ### read header 
    sampleLine = fhs.readline()
    sampleLine = sampleLine.rstrip()
    sampleList = sampleLine.split(sep="\t")
    tt = sampleList.pop(0)
    tt = sampleList.pop(0)

    indexList = []
    foundPlateDict={}
    sampleCount = len(sampleList)
    plateWell2Index = {}
    
    for ii in range(sampleCount):
        sampleName = sampleList[ii]
        if (mode=="generic_sample"):
            if (sampleName in sampleNameList) and (sampleName not in plateWell2Index):
                plateWell2Index[sampleName] = ii
        else:            
            if ("__" not in sampleName):
                continue
            plateWell = sampleName.split("__")[1]
            m= re.match("(.+)_\w{3}$", plateWell)
            if m:
                plate = m[1]
            else:
                print(f"Warning: sample platewell {plateWell} is skipped")
                continue    
            
            if (mode=="sample") and (plateWell in plateWellList) and (plateWell not in plateWell2Index):
                plateWell2Index[plateWell] = ii
                
            if (mode=="plate") and (plate in plateList) and (plateWell not in plateWell2Index):
                indexList.append(ii)
                plateWell2Index[plateWell] = ii
                if plate in foundPlateDict:
                    foundPlateDict[plate]+=1
                else:
                    foundPlateDict[plate]=1

        
        
    #first make sure the sample list or plate list provide in familyFile are  present in the hap_genotype file
    if (mode=="generic_sample"):
        missingSamples = []
        for t in sampleNameList:          
            if t in plateWell2Index:
                indexList.append(plateWell2Index[t])
            else:
                missingSamples.append(t)
        if len(missingSamples) >0:
            print (f"The following individuals are not found in hap_genotype file: Please correct them and try again.")
            print ("\n".join(missingSamples))
            sys.exit()
            
    if (mode=="sample"):
        missingSamples = []
        for t in plateWellList:          
            if t in plateWell2Index:
                indexList.append(plateWell2Index[t])
            else:
                missingSamples.append(t)
        if len(missingSamples) >0:
            print (f"The following individuals are not found in hap_genotype file: Please correct them and try again.")
            print ("\n".join(missingSamples))
            sys.exit()

    
    if (mode=="plate"):
        missingPlateList=[]
        print ("Found samples per plate in hap_genotype:")
        for p in plateList:
            if p in foundPlateDict:
                print(f"{p}\t{foundPlateDict[p]}")                 
            else:
                missingPlateList.append(p)
                print(f"{p}\t0")                   
        if len(missingPlateList) >0:
            print (f"The following plates are not found in hap_genotype file: Please correct them and try again.")
            print ("\n".join(missingPlateList))
            sys.exit()

    filePath =  f"{args.output}/hap_genotype"
    fh = open (filePath,'w')
    fh.write ("Locus\tHaplotypes")
    
    for sIndex in indexList:        
        sampleName = sampleList[sIndex]
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
            pass
            #print(f"{markerName}: no data")

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
            pass
            #print(f"{markerName}: no used alleles")

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


