#!/usr/bin/env python3
import sys
import argparse
import os
from os.path import isfile, join
import re

import pandas as pd
from functools import reduce
import glob

def read_tab_delimited_files(file_paths):
    tables = []
    for file_path in file_paths:
        df = pd.read_csv(file_path, delimiter='\t', header=0)
        df = df.drop(columns=['Haplotypes'], errors='ignore')
        tables.append(df)
        sampleCount =len(df.columns)-1
        print(f"{sampleCount} samples in {file_path}")
    return tables

def merge_tables(tables):
    merged_table = reduce(lambda left, right: pd.merge(left, right, on="Locus", how='outer'), tables)
    merged_table.fillna('./.:0', inplace=True)
    return merged_table

    
    

def calc_allele_frq(row):
    totalGametes =0
    alleleCounts = {}
    newHapStr = "";
    for gtstr in row[1:]:
        GT = gtstr.split(":")                
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
    if (totalGametes>0):
        sortedAlleles = sorted(alleleCounts.items(), key=lambda kv: kv[1], reverse=True)
        # re-calculate allele frq based on used allelels
        for alleleCount in sortedAlleles:
            aId = alleleCount[0]
            b = "%4.2f" %(int(alleleCount[1])/totalGametes)
            newHapStr += f"{aId}({b});"  
    return (newHapStr)
        
def keep_first_columns(df):
    # Identify the columns to keep
    cols_to_keep = ~df.columns.duplicated(keep='first')
    # Filter the DataFrame to keep only these columns
    return df.loc[:, cols_to_keep]   
        

def main():
    if sys.version_info[0] < 3:
        raise Exception("This code requires Python 3.")

    parser = argparse.ArgumentParser(description='Run slicer.')
    parser.add_argument('-i','--input',type=str,required=True,help='Input files, it must be the output hap_genotype files from the amplicon.py tool, separated by comma')
    parser.add_argument('-o','--output',type=str,required=True,help='Output hap_genotype file.')


    args=parser.parse_args()



    #process hap_genotype file
    input_file_list = args.input.split(",")
    if (len(input_file_list) <2):
        print (f"Error: There must be two or more input files, separated by comma. Your input file string is {args.input}")
        sys.exit()

    for file in input_file_list:
        if (not os.path.isfile(file)):
            print(f"Error: the input file {file} does not exist!")
            sys.exit()

    myTables = read_tab_delimited_files(input_file_list)
    myMerged = merge_tables(myTables)
    myMerged.insert(1, 'Haplotypes', myMerged.apply(calc_allele_frq, axis=1))
    
    #cols_to_keep = ~myMerged.columns.duplicated(keep='first')
    # Filter the DataFrame to keep only these columns
    #myMerged = myMerged.loc[:, cols_to_keep]   
    duplicated_cols = myMerged.columns[myMerged.columns.duplicated(keep='first')]
    if (len(duplicated_cols)>0):
        print ("These columns will be dropped")
        print (duplicated_cols)
    cols_to_keep = ~myMerged.columns.duplicated(keep='first')
    myMerged = myMerged.loc[:, cols_to_keep]
    sampleCount =len(df.columns)-2
    print (f"{sampleCount} samples in merged table. ")
    myMerged.to_csv(args.output, sep='\t', index=False)



if __name__=="__main__":
    main()


