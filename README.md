# amplicon.sh

This python script (amplicon.py) is developed for the Vitisgen2 project (https://www.vitisgen2.org/) to analyzes multi-plexed amplicon sequencing data, optimized for high-throughput IDT rhAmpSeq data. Another script (to_lep_map.pl) converts the output from amplicon.py to a VCF file that can be loaded into LepMap3 software for genetic linkage map construction. The alleles in the VCF file (A, C, G, T) are not actual nucleotide alleles, but symbols representing up to 4 haplotytpe alleles per marker. There is a lookup table in the output files to give you the actual sequences of each haplotype represented by the "ACGT".

## Getting Started


### Prerequisites
1. Python 3 (tested on python 3.6) and PERL (tested on 5.22)
2. Python module: BioPython https://biopython.org/
3. Other software. 
The following commands should be installed and in the PATH:
* bbmerge.sh: https://sourceforge.net/projects/bbmap/
* cutadapt: https://cutadapt.readthedocs.io/en/stable/
* muscle (optional): https://www.drive5.com/muscle/

### Installation
Download the two scripts: amplicon.py and to_lep_map.pl, and put them in any directory

### Usage
1. Preparing data files.
  Put the following items in the same directory:
  * Sample file: A tab-delimited text file with three columns. 1)Sample Name; 2) Paired-end sequence file 1 (fastq or fastq.gz); 3) Paired-end sequence file 2.
  * Key file: A tab delimited text file with three columns. 1) marker name; 2) 5' PCR primer sequence; 3) 3' PCR primer sequence.
  * All fastq.gz files listed in the sample file.

2. In the data directory, execute this command:
amplicon.sh -s sampleFileName -k keyFileName -o outputDirName -j 10 -a 0.15
 * -j 10:  process 10 samples simultaneously. This should not exceed the total number of CPU cores on your computer, and should not exceed 20 even if you have more cores du to IO limitation.
 * =a 0.15: minimum minor allele frequency

### Other parameters
  *-h	show this help message and exit
  *-i	Skip steps. e.g. "-i 12" to skip steps 1 and 2. the steps are: 1. split reads by primers; 2. identify haplotypes across population, and optionally run PCR error correction if set "-e 1"; 3 call genotypes
  *-j	Number of simultaneous jobs. Default:8
  *-t	Number of threads per job. Default:4
  *-c	Minimum number of samples per haplotypes. Default:10
  *-n	Maximum number of unique haplotypes per sample to be kept in the first step. Default:20
  *-m	Maximum number of haplotypes per marker in the population. Default:1000
  *-a	Minimum minor allele frequency Default:0.05
  *-l	Minimum haplotype length (after removing the primers. It must be an integer 1 or larger.) Default:20
  *-d	Whether to merge the duplicate samples. 1: Merge; 0: Do not merge and the duplicated sample will be named <sampleName>__<index starting from 1>. Default:1
  -e	Correct PCR errors based on allele frequency (only applicable for biparental families). 0: No correction; 1: Correct error in bi-parental population based on allele read count distribution in the population. Default:0, no correction
  -p	Ploidy, default 2


## Authors
* **Qi Sun**
* **Cheng Zou**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments
