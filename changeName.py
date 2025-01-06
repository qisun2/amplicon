#!/usr/bin/env python3
import pandas as pd
import sys
import argparse
import os
def main():
    if sys.version_info[0] < 3:
        raise Exception("This code requires Python 3.")

    parser = argparse.ArgumentParser(description='Run slicer.')
    parser.add_argument('-i','--input',type=str,required=True,help='Input hap_genotype file name.')
    parser.add_argument('-o','--output',type=str,required=True,help='Output hap_genotype file name.')
    parser.add_argument('-n','--nameFile',type=str,required=True,help='sample name file. It must be a tab-delimited text file, with two columns, no header. The first column is the original name, the second column is the new name.')
    

    args=parser.parse_args()

        
    if (args.input != None) and (not os.path.isfile(args.input)):
        parser.print_usage()
        print(f"Error: input file {args.input} does not exist!")
        sys.exit()

    if (args.nameFile != None) and (not os.path.isfile(args.nameFile)):
        parser.print_usage()
        print(f"Error: nameFile file {args.nameFile} does not exist!")
        sys.exit()

    df = pd.read_csv(args.input, delimiter='\t', dtype=str, header=0)
    column_mapping = pd.read_csv(args.nameFile, delimiter='\t', dtype=str, header=None, names=['old_name', 'new_name'])
    
    # Create a dictionary from the second table to map old names to new names
    column_rename_dict = dict(zip(column_mapping['old_name'], column_mapping['new_name']))

    # Replace the column names in the first DataFrame
    df.rename(columns=column_rename_dict, inplace=True)


    # Save the updated DataFrame to a new file if needed
    df.to_csv(args.output, sep='\t', index=False)



if __name__=="__main__":
    main()


