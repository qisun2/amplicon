# amplicon.py

This python script (amplicon.py) is developed for the Vitisgen2 project (https://www.vitisgen2.org/) to analyzes multi-plexed amplicon sequencing data, optimized for high-throughput IDT rhAmpSeq data. Another script (to_lep_map.pl) converts the output from amplicon.py to a VCF file that can be loaded into LepMap3 software for genetic linkage map construction. The alleles in the VCF file (A, C, G, T) are not actual nucleotide alleles, but symbols representing up to 4 haplotytpe alleles per marker. There is a lookup table in the output files to give you the actual sequences of each haplotype represented by the "ACGT".  

If you want to process your amplicon sequencing data using GATK to call variants, you can use this code: https://bitbucket.org/cornell_bioinformatics/gatk4amplicon  

## Getting Started


### Prerequisites
1. Python 3 (tested on python 3.6) and PERL (tested on 5.22)
2. Python module: BioPython https://biopython.org/
3. Other software. 
The following commands should be installed and in the PATH:
* bbmerge.sh: https://sourceforge.net/projects/bbmap/
* cutadapt: https://cutadapt.readthedocs.io/en/stable/ (v3 or above)
* muscle (optional): https://www.drive5.com/muscle/
* python package swalign (required to run amplicon_dada2.py , https://pypi.org/project/swalign/ )
* R package dada2 (required to run amplicon_dada2.py , https://bioconductor.org/packages/release/bioc/html/dada2.html )

### Installation
Download the two scripts: amplicon.py and to_lep_map.pl, and put them in any directory

### Usage
1. Preparing data files.  
    Put the following items in the same directory:  
  * Sample file: A tab-delimited text file with four columns. 1)Sample Name; 2)sample plate_well, it can be any string to uniquely identify a sample if the sample names are duplicated; 3) Paired-end sequence file 1 (fastq or fastq.gz); 4) Paired-end sequence file 2.  
  * Key file: A tab delimited text file with three columns. 1) marker name; 2) 5' PCR primer sequence; 3) 3' PCR primer sequence.  
  * All fastq.gz files listed in the sample file.  

2. In the data directory, execute this command:  
amplicon.py -s sampleFileName -k keyFileName -o outputDirName -j 10 -a 0.15  
 * -j 10:  process 10 samples simultaneously. This should not exceed the total number of CPU cores on your computer, and should not exceed 20 even if you have more cores du to IO limitation.  
 * -a 0.15: minimum minor allele frequency  

### Output files
  * hap_genotype: A matrix with all genotypes. Each row is a marker, each column is a sample.
  * hap_genotype_matrix: A file transposed from hap_genotype file, and without read count per allele information. 
  * HaplotypeAllele.fasta: A fasta file with haplotype sequence per allele.
  * topHaplotypeAllele.fasta: A fasta file with the top allele per marker


### Other parameters
  * -h	show this help message and exit
  * -i	Skip steps. e.g. "-i 12" to skip steps 1 and 2. the steps are: 1. split reads by primers; 2. identify haplotypes across population, and optionally run PCR error correction if set "-e 1"; 3 call genotypes
  * -j	Number of simultaneous jobs. Default:8
  * -t	Number of threads per job. Default:1
  * -c	Minimum number of samples per haplotypes. Default:10
  * -n	Maximum number of unique haplotypes per sample to be kept in the first step. Default:20
  * -m	Maximum number of haplotypes per marker in the population. Default:1000
  * -a	Minimum minor allele frequency Default:0.01
  * -l	Minimum haplotype length (after removing the primers. It must be an integer 1 or larger.) Default:50
  * -d	Whether to merge the duplicate samples. 1: Merge; 0: Do not merge and the duplicated sample will be named <sampleName>__<index starting from 1>. Default:1
  * -e	Correct PCR errors based on allele frequency (only applicable for biparental families). 0: No correction; 1: Correct error in bi-parental population based on allele read count distribution in the population. Default:0, no correction
  * -p	Ploidy, default 2
  * -r	Maximum read count ratio between the two alleles in each sample, default 20
  * -z	Mismatch rate between pcr primer and reads, default 0.1

### to_lep_map.pl script
As many software, e.g. Lep-MAP3, requires vcf file format. This script is provided to convert the hap_genotype from previous step to a "fake" VCF file, in which up to 4 haplotype alleles are represented with A C G T codind.  There is a lookup table to show the corresponding haplotype allele sequences for "A" "C" "G" and "T".  
  *  Usage:  to_lep_map.pl -g hap_genotype -f minorAlleleFrequency -b blankSample -m maternalSample -p paternalSample -l marker2pos -n familyName  
  *  -b, -m, -p, -l are optional. "-b" "-m" and "-p" are integer index of the blank, maternal and paternal samples in your sample list. If multiple samples, separate the index with comma. The index should be 1-based, so that the first sample is 1. "-l" is to specify the physical positions of the markers. It should be a tab-delimited table with 3 columns: markerName, chr, pos. 

### hap2realVCF tool
This tool can convert the amplicon.py output into a real VCF file (not like the fake vcf file from to_lep_map.pl script). The SNPs/indels in the vcf file are not phased in haplotypes. See README.md file in the directory how to setup and run this script.

### merge.py script
To merge multiple hap_genotype files, run the command: merge.py -i inputfile1,inputfile2,inputfile3 -o outputfile 
(The input files must be the hap_genotype files from amplicon.py. If there are duplcated sample names, only firt ones are kept)

### slice.py script
You can slice by either a list of plates, or a list of individuals.

To slice by a list of individuals, prepare a sample name file, with one sample name per line. The sample name should be in the format: plateName_wellName, e.g. vDNAcad794B07_E03. Any text after first column (tab-delimited) will be ignored.
Run the command: slice.py -i inputDirectory -o outputDirectory -f sampleFileName -a MAF

(In the output file, the samples will be ordered based on what is in the file)

To slice out a family by a list of plates, prepare a plate name file, with one plate name per line.
Run the command: slice.py -i inputDirectory -o outputDirectory -p plateFileName -a MAF 

This new code will slice out samples from the genotyping matrix. The inputDirectory should contain hap_genotype and optionally markerToSampleReadCountMatrix from the amplicon.py or merge.py script. MAF is the minimum allele frequency in the family). After run the code, you should see a new output directory, with hap_genotype and markerToSampleReadCountMatrix (if provided as input). The hap_genotype file can be used as input for to_lep_map.pl to convert to VCF files and lepmap pedigree files.

### changeName.py script
changeName.py -i input_hap_genotype_filename -o output_hap_genotype_filename -n namefile

The name file is a tab delimited text file with two columns (no header line). oldname and newname.

## Authors
* **Qi Sun**
* **Cheng Zou**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments