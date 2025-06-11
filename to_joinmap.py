#!/usr/bin/env python3
import pandas as pd
import sys
import argparse
import os
import re


def main():
    if sys.version_info[0] < 3:
        raise Exception("This code requires Python 3.")

    parser = argparse.ArgumentParser(description='Run slicer.')
    parser.add_argument('-i','--input',type=str,required=True,help='Input hap_genotype file name.')
    parser.add_argument('-o','--output',type=str,required=True,help='Output hap_genotype file name.')
    parser.add_argument('-n', '--familyname', required=False, type=str, default='mypopulation', help='Name of the family')
    parser.add_argument('-m','--maternal',type=str,required=True,help='Maternal parents.')
    parser.add_argument('-p','--paternal',type=str,required=True,help='Paternal parents.')
    parser.add_argument('-a','--abxcc',type=str, default="1", required=False,help='Options to code abcc markers. 1:abbb; 2:abcd; 3:abbb and aaab')
    parser.add_argument('-x','--missing',type=float,required=False, default=0.5, help='Missing rate.')
    global args
    args=parser.parse_args()

        
    if (args.input != None) and (not os.path.isfile(args.input)):
        parser.print_usage()
        print(f"Error: input file {args.input} does not exist!")
        sys.exit()

    
    maternal_list = int_list = [int(x)+1 for x in re.findall(r'\d+', args.maternal)]
    paternal_list = int_list = [int(x)+1 for x in re.findall(r'\d+', args.paternal)]
    
    
    if (len(maternal_list)==0):
        print(f"The maternal list must be a list of integers")
        sys.exit()
    if (len(paternal_list)==0):
        print(f"The paternal list must be a list of integers")
        sys.exit()

    df = pd.read_csv(args.input, delimiter='\t', dtype=str, header=0)
    



    maternal_names = df.columns[maternal_list].tolist()
    paternal_names = df.columns[paternal_list].tolist()
    
    #progenies
    all_indices = range(2, len(df.columns))
    exclude_set = set(maternal_list + paternal_list)
    progeny_list = [x for x in all_indices if x not in exclude_set]
    progeny_names_ori= df.columns[progeny_list].tolist()
    
    progeny_names=[]
    for pname in progeny_names_ori:
        pname, plate = pname.split("__")
        if len(pname)>20:
            pname=pname[:20]
        progeny_names.append(pname)
        
    WHLOG = open (args.output + ".run.log", "wt")
    WHINFO = open (args.output + ".marker_info", "wt")
    
    WHLOG.write (f"Mathernal parents ({len(maternal_names)}):\n\t{maternal_names}\n")
    WHLOG.write (f"Paternal parents ({len (paternal_names)}):\n\t{paternal_names}\n")
    WHLOG.write (f"Progenies ({len(progeny_list)}):\n\t{progeny_names}\n")
    
    coding_stat = {}
    marker_stat = {"PASSED":0}
    output_gt_list = []
    WHINFO.write(f"Marker\tOriCode\tMaternalAllele\tPaternalAllele\tMaternalGTs\tPaternalGTs\tPASS\tFilterCode\tCodeInLoc\t#MissingOneAllele\t#MissingTwoAllele\t#MissingRate\t" + "\t".join(progeny_names) + "\n")
    WHLOG.write (f"\n\nMarkers filtered due to high missing rate (cutoff: {args.missing})\nMarker\tCode\t#Missing One allele\tMissing Two alleles\tOverall Missing rate\n")
    for index, row in df.iterrows():
        maternal_gts = [row[i] for i in df.columns[maternal_list]]
        paternal_gts = [row[i] for i in df.columns[paternal_list]]
        progeny_gts = [row[i] for i in df.columns[progeny_list]]
        
        markerName = row.iloc[0]
        markerInfo = row.iloc[1]
        
        
        markerName = markerName[:20] if len(markerName) > 20 else markerName
        
        maternal_alleles = get_top_two_alleles(maternal_gts)
        #print (maternal_alleles, maternal_gts)
        paternal_alleles = get_top_two_alleles(paternal_gts)
        #print (paternal_alleles, paternal_gts)
        
        code = "unknown"
        
        if (len(maternal_alleles) == 0) or (len(paternal_alleles) == 0):
            code = "Parental_gt_unk"
        elif (len(maternal_alleles) == 1) and (len(paternal_alleles) == 1):
            code= "aabb"
        elif (len(maternal_alleles) == 2) and (len(paternal_alleles) == 2):
            if (maternal_alleles[0] in paternal_alleles) and (maternal_alleles[1] in paternal_alleles):
                code = "hkxhk"
            elif (maternal_alleles[0] not in paternal_alleles) and (maternal_alleles[1] not in paternal_alleles):
                code = "abxcd"
            else:
                if (maternal_alleles[0] not in paternal_alleles):
                    t = maternal_alleles[0]
                    maternal_alleles[0]= maternal_alleles[1]
                    maternal_alleles[1] =t
                if (paternal_alleles[0] not in maternal_alleles):
                    t = paternal_alleles[0]
                    paternal_alleles[0]= paternal_alleles[1]
                    paternal_alleles[1] =t                    
                code = "efxeg"
        elif (len(maternal_alleles) == 2) and (len(paternal_alleles) == 1):
            if (paternal_alleles[0] in maternal_alleles):
                if paternal_alleles[0] == maternal_alleles[1]:
                    t = maternal_alleles[0]
                    maternal_alleles[0]= maternal_alleles[1]
                    maternal_alleles[1] =t                    
                code = "lmxll"
            else:
                code = "abxcc"
        elif (len(maternal_alleles) == 1) and (len(paternal_alleles) == 2):
            if maternal_alleles[0] in paternal_alleles:               
                if maternal_alleles[0] == paternal_alleles[1]:
                    t = paternal_alleles[0]
                    paternal_alleles[0]= paternal_alleles[1]
                    paternal_alleles[1] =t
                code = "nnxnp"
            else:
                code = "ccxab"        
        
        if code in coding_stat:
            coding_stat[code] += 1
        else:
            coding_stat[code] =1
        
        ExtraMarkerInfoStr = ""
        
        if code in ("abxcd", "hkxhk", "efxeg", "lmxll", "nnxnp", "abxcc", "ccxab"):
            code_in_locfile = code
            
            gtline_extra = ""
            
            if code in ("abxcc", "ccxab"):
                if (args.abxcc=="1") or (args.abxcc=="3"):
                    if code == "abxcc":
                        code_in_locfile = "lmxll"
                    elif code=="ccxab":
                        code_in_locfile = "nnxnp"
                elif args.abxcc=="2":
                    code_in_locfile = "abxcd"
                    
                if args.abxcc == "3" and code in ("abxcc", "ccxab"):
                    markerName = markerName[:18]
                    m1 = markerName + "-a"
                    m2 = markerName + "-b"
                    gtline = f"{m1.ljust(25, ' ')}  <{code_in_locfile}>\n "
                    gtline_extra = f"{m2.ljust(25, ' ')}  <{code_in_locfile}>\n "
                else:
                    gtline = f"{markerName.ljust(25, ' ')}  <{code_in_locfile}>\n "
            else:
                gtline = f"{markerName.ljust(25, ' ')}  <{code_in_locfile}>\n "
               
            evline = ""            
            missingOneAlleleCount = 0
            missingTwoAlleleCount = 0
            total = len (progeny_gts)
            

            for progeny_gt in progeny_gts:            
                gt, evidence = convert_gt(progeny_gt, maternal_alleles, paternal_alleles, code)
                if gt == "--":
                    missingTwoAlleleCount +=1
                elif "-" in gt:
                    missingOneAlleleCount += 1
               
                if (code in ("abxcc", "ccxab")) and (args.abxcc == "3"):  
                    gt1, gt2 = gt.split("&&")
                    gtline = gtline + " " + gt1
                    gtline_extra = gtline_extra + " " + gt2
                else:
                    gtline = gtline + " " + gt
                    
                evline = evline + evidence + "\t"
                
            missingRate = (missingOneAlleleCount + missingTwoAlleleCount)/total
            ExtraMarkerInfoStr += f"{missingOneAlleleCount}\t{missingTwoAlleleCount}\t{missingRate}\t{evline}"
            
            if missingRate < args.missing:                
                if code in ("ccxab", "abxcc"):
                    extra_info = f"  Note: SEG type converted from {code}"
                    output_gt_list.append(gtline + f" ; maternal alleles: {maternal_alleles} paternal alleles: {paternal_alleles} {extra_info}" )
                    if args.abxcc == "3":
                        output_gt_list.append(gtline_extra + f" ; maternal alleles: {maternal_alleles} paternal alleles: {paternal_alleles} {extra_info}" )
                else:
                    output_gt_list.append(gtline + f" ; maternal alleles: {maternal_alleles} paternal alleles: {paternal_alleles}" )
                ExtraMarkerInfoStr=f"Passed\t\t{code_in_locfile}\t" + ExtraMarkerInfoStr
                  
                marker_stat["PASSED"] += 1
            else:
                ExtraMarkerInfoStr=f"Removed\tMissingRate\t{code_in_locfile}\t" + ExtraMarkerInfoStr
                WHLOG.write(f"{markerName}\t{code}\t{missingOneAlleleCount}\t{missingTwoAlleleCount}\t{missingRate}\n")
                
                ccc = "NOT PASSED missing rate filter " + code
                if ccc in marker_stat:
                    marker_stat[ccc] += 1
                else:
                    marker_stat[ccc] = 1
        else:
            filterReason = "CodeNotUsed"
            if code == "Parental_gt_unk":
                filterReason = "NoParentalAlleles"            
            ExtraMarkerInfoStr=f"Removed\t{filterReason}\t{code_in_locfile}\t"
        WHINFO.write(f"{markerName}\t{code}\t{', '.join(maternal_alleles)}\t{', '.join(paternal_alleles)}\t{'--'.join(maternal_gts)}\t{'--'.join(paternal_gts)}\t{ExtraMarkerInfoStr}\n")
    
    WHLOG.write("\n\nCode stats:\n")
    for code in coding_stat:
        WHLOG.write(f"{code}\t{coding_stat[code]}\n")


    WHLOG.write("\n\nMarker stats:\n")
    for s in marker_stat: 
        WHLOG.write(f"{s}\t{marker_stat[s]}\n")    

        
    WHLOG.close()
    WHINFO.close()

    
    #write to output file
    gt_output = "\n".join(output_gt_list)
    indiv_names = "\n".join(progeny_names)
    
    
    text_block = f"""; converted from amplicon genotyping matrix
    
name = {args.familyname}
popt = CP
nloc = {len(output_gt_list)}
nind = {len(progeny_list)}

{gt_output}

individual names:
{indiv_names}

"""


    WHOUT = open (args.output + ".loc", "wt")
    WHOUT.write(text_block)
    WHOUT.close()
    
def convert_gt (progeny, mother, father, code):
    alleles, readcounts = progeny.split(":")
    allele_list = alleles.split("/")
    readcount_list = readcounts.split(",")
    
    if code == "abxcd":
        g1 = "-"
        g2 = "-"
        if (mother[0] in allele_list):
            g1="a"
        elif (mother[1] in allele_list):
            g1="b"
        if (father[0] in allele_list):
            g2="c"
        elif (father[1] in allele_list):
            g2="d"
            
        return [f"{g1}{g2}", f"{alleles}={g1}{g2}"]
    elif code == "hkxhk":
        g1 = "-"
        g2 = "-"
        if (mother[0] in allele_list):
            g1="h"
        if (mother[1] in allele_list):
            g2="k"
        if g1 == "h" and g2 =="-":
            g2 = "h"
        if g1 == "-" and g2 =="k":
            g1 = "k"
        return [f"{g1}{g2}", f"{alleles}={g1}{g2}"]
    elif code == "efxeg":
        gt = "--"
        if (mother[0] in allele_list) and (father[1] in allele_list):
            gt = "eg"
        elif (mother[0] in allele_list) and (mother[1] in allele_list):
            gt = "ef"
        elif (mother[1] in allele_list) and (father[1] in allele_list):
            gt = "fg"
        elif (mother[0] in allele_list):
            gt = "ee"
        return [f"{gt}", f"{alleles}={gt}"]
    elif code == "lmxll":
        gt = "--"
        if (mother[0] in allele_list) and (mother[1] in allele_list):
            gt = "lm"
        elif (mother[0] in allele_list) and (mother[1] not in allele_list):
            gt = "ll"
        return [f"{gt}", f"{alleles}={gt}"]
    elif code == "nnxnp":
        gt = "--"
        if (father[0] in allele_list) and (father[1] in allele_list):
            gt = "np"
        elif (father[0] in allele_list) and (father[1] not in allele_list):
            gt = "nn"
        return [f"{gt}", f"{alleles}={gt}"]
    elif code == "abxcc":
        if args.abxcc == "1":
            gt = "--"        
            if (mother[0] not in allele_list) and (mother[1] not in allele_list):
                return ["--", f"{alleles}=--"]
                
            if mother[0] in allele_list:
                gt = "lm"
            else:
                gt = "ll"
            return [gt, f"{alleles}={gt}"]
        elif args.abxcc == "2":
            gt = "--"        
            if (mother[0] not in allele_list) and (mother[1] not in allele_list):
                return ["--", f"{alleles}=--"]
                
            if (mother[0] in allele_list) and (mother[1] not in allele_list):
                gt = "ac"
            elif (mother[1] in allele_list) and (mother[0] not in allele_list):
                gt = "bc"
            else:
                gt = "ab"
            return [gt, f"{alleles}={gt}"]
        elif args.abxcc == "3":
            gt1="--"
            gt2="--"        
            if (mother[0] not in allele_list) and (mother[1] not in allele_list):
                return ["--&&--", f"{alleles}=--"]
            

            if (mother[0] in allele_list):
                gt1 = "lm"
            else:
                gt1 = "ll"
                
            if (mother[1] in allele_list):
                gt2 = "lm"
            else:
                gt2 = "ll"
                
            return [f"{gt1}&&{gt2}", f"{alleles}={gt1}&&{gt2}"]
    elif code == "ccxab":
        if args.abxcc == "1":
            gt = "--"        
            if (father[0] not in allele_list) and (father[1] not in allele_list):
                return ["--", f"{alleles}=--"]
                
            if father[0] in allele_list:
                gt = "lm"
            else:
                gt = "ll"
            return [gt, f"{alleles}={gt}"]
        elif args.abxcc == "2":
            gt = "--"        
            if (father[0] not in allele_list) and (father[1] not in allele_list):
                return ["--", f"{alleles}=--"]
                
            if (father[0] in allele_list) and (father[1] not in allele_list):
                gt = "ca"
            elif (father[1] in allele_list) and (father[0] not in allele_list):
                gt = "cb"
            else:
                gt = "ab"
            return [gt, f"{alleles}={gt}"]
        elif args.abxcc == "3":
            gt1="--"
            gt2="--"        
            if (father[0] not in allele_list) and (father[1] not in allele_list):
                return ["--&&--", f"{alleles}=--"]
            

            if (father[0] in allele_list):
                gt1 = "lm"
            else:
                gt1 = "ll"
                
            if (father[1] in allele_list):
                gt2 = "lm"
            else:
                gt2 = "ll"
                
            return [f"{gt1}&&{gt2}", f"{alleles}={gt1}&&{gt2}"]
    else:
        return ["--", alleles]
            

def get_top_two_alleles (gt_list):
    alleleToReadCounts = {}
    for gt_str in gt_list:
        alleles, readcounts = gt_str.split(":")
        allele_list = alleles.split("/")
        readcount_list = readcounts.split(",")
        
        if allele_list[0] == ".":
            continue
            
        if allele_list[0] == allele_list[1]:
            allele_list.pop()
            
        for a in allele_list:
            readcount = readcount_list.pop(0)
            if a in alleleToReadCounts:
                alleleToReadCounts[a] += int(readcount)
            else:
                alleleToReadCounts[a] =int(readcount)
    top2alleles = sorted(alleleToReadCounts, key=alleleToReadCounts.get, reverse=True)[:2]
    if (len(top2alleles) >1) and ((alleleToReadCounts[top2alleles[0]]/alleleToReadCounts[top2alleles[1]]) > 10):
        top2alleles.pop()
    return (top2alleles)
        
            
        


if __name__=="__main__":
    main()


