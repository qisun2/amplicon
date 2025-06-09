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
    parser.add_argument('-x','--missing',type=float,required=False, default=0.5, help='Missing rate.')
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
    progeny_names= df.columns[progeny_list].tolist()
    
    WHLOG = open (args.output + ".log", "wt")
    WHINFO = open (args.output + ".marker_info", "wt")
    
    WHLOG.write (f"Mathernal parents: {maternal_names}\n")
    WHLOG.write (f"Paternal parents: {paternal_names}\n")
    WHLOG.write (f"Progenies {len(progeny_list)}: {progeny_names}\n")
    
    WHTEST = open (args.output + "testing.log", "wt")
    
    output_gt_list = []
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
            code = "missing"
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

        WHINFO.write(f"{markerName}\t{code}\t{'/'.join(maternal_alleles)}\t{'/'.join(paternal_alleles)}\t{'--'.join(maternal_gts)}\t{'--'.join(paternal_gts)}\n")
        if code in ("abxcd", "hkxhk", "efxeg", "lmxll", "nnxnp", "abxcc", "ccxab"):
            
            if code == "abxcc":
                if len(markerName)>18:
                    markerName = markerName[:20]
                m1 = markerName + "-1"
                m2 = markerName + "-2"
                
                gtline = f"{m1.ljust(25, ' ')}  <lmxll>\n "
                gtline_extra = f"{m2.ljust(25, ' ')}  <lmxll>\n "           
            elif code=="ccxab":
                if len(markerName)>18:
                    markerName = markerName[:20]
                m1 = markerName + "-1"
                m2 = markerName + "-2"
                gtline = f"{m1.ljust(25, ' ')}  <nnxnp>\n "
                gtline_extra = f"{m2.ljust(25, ' ')}  <nnxnp>\n "               
            else:
                gtline = f"{markerName.ljust(25, ' ')}  <{code}>\n "
                
            evline = f"{code}\t{markerName}\t{maternal_alleles}\t{paternal_alleles}"            
            missingCount = 0
            total = len (progeny_gts)
            

            for progeny_gt in progeny_gts:            
                gt, evidence = convert_gt(progeny_gt, maternal_alleles, paternal_alleles, code)
                if "-" in gt:
                    missingCount = missingCount + 1
                    
                if code in ("ccxab", "abxcc"):
                    gt1, gt2 = gt.split("&&")
                    gtline = gtline + " " + gt1
                    gtline_extra = gtline_extra + " " + gt2                 
                else:
                    gtline = gtline + " " + gt 
                
                evline = evline + "\t" + evidence
            if (missingCount/total)< args.missing:
                if code in ("ccxab", "abxcc"):
                    output_gt_list.append(gtline + f" ; maternal alleles: {maternal_alleles} paternal alleles: {paternal_alleles} {code}-1" )
                    output_gt_list.append(gtline_extra + f" ; maternal alleles: {maternal_alleles} paternal alleles: {paternal_alleles} {code}-2" )
                else:
                    output_gt_list.append(gtline + f" ; maternal alleles: {maternal_alleles} paternal alleles: {paternal_alleles}" )
                WHTEST.write(f"{evline}\n")
            else:
                WHLOG.write(f"missing\t{markerName}\t{missingCount/total}\n")
    WHLOG.close()
    WHTEST.close()

    #write to output file
    gt_output = "\n".join(output_gt_list)
    
    indiv_names = ""
    for pname in progeny_names:
        pname, plate = pname.split("__")
        if len(pname)>20:
            pname=pname[:20]
        indiv_names = indiv_names + pname + "\n"
    
    
    text_block = f"""; converted from amplicon genotyping matrix
    
name = {args.familyname}
popt = CP
nloc = {len(output_gt_list)}
nind = {len(progeny_list)}

{gt_output}

individual names:
{indiv_names}
"""


    WHOUT = open (args.output, "wt")
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
        gt_a = "--"
        gt_b = "--"
        
        if (mother[0] not in allele_list) and (mother[1] not in allele_list):
            return ["--&&--", f"{alleles}=--"]
            
        if mother[0] in allele_list:
            gt_a = "lm"
        else:
            gt_a = "ll"
            
        
        if mother[1] in allele_list:
            gt_b = "lm"
        else:
            gt_b = "ll"
            
        return [f"{gt_a}&&{gt_b}", f"{alleles}={gt_a}&&{gt_b}"]
    elif code == "ccxab":
        gt_a = "--"
        gt_b = "--"
        
        if (father[0] not in allele_list) and (father[1] not in allele_list):
            return ["--&&--", f"{alleles}=--"]
            
        if father[0] in allele_list:
            gt_a = "np"
        else:
            gt_a = "nn"
            
        
        if father[1] in allele_list:
            gt_b = "np"
        else:
            gt_b = "nn"
            
        return [f"{gt_a}&&{gt_b}", f"{alleles}={gt_a}&&{gt_b}"]
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


