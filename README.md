This packge is developed for the Vitisgen2 & Vitisgen3 projects ([https://vitisgen3.umn.edu](https://vitisgen3.umn.edu)) to analyze multiplexed amplicon sequencing data. The output is a multi-allelic haplotype genotyping matrix. If a known-haplotype allele sequence file is provided, the script can be used for genotyping with known alleles. 

### Prerequisites

1. Python 3 (tested on python 3.6 and 3.9) and PERL (tested on 5.22 and 5.32)
2. Python module: BioPython https://biopython.org/
3. Other software. 
   The following commands should be installed and in the $PATH:

* cutadapt (required): https://cutadapt.readthedocs.io/en/stable/ (v3 or above)
* bbmap (required for paired-end reads only): https://sourceforge.net/projects/bbmap/


### Installation

Download the  repository: 

```
git clone https://bitbucket.org/cornell_bioinformatics/amplicon.git
```



### Applications

#### 1. Call genotypes from raw sequencing data 

**Script:** amplicon.py

**Steps:**

1. **<u>Preparing input data files.</u>**  

  * **Sample file:** A tab-delimited text file with either three (single-end reads) or four columns (paired-end reads). No header line.
    * Sample Name;
    * PlateName_Well of the sample (eg. vDNAcad794B07_E03. If not available, duplicate the sample names from the first column);
    * Path of paired-end sequence file 1 (fastq or fastq.gz);
    * Path of paired-end sequence file 2 (skip for single end data) .  

  * **Key file:** A tab delimited text file with three columns. No header line. 1) marker name; 2) 5' PCR primer sequence; 3) 3' PCR primer sequence. 
    * By default, the sequences of the two primers should be in F-R orientation. If the orientation is F-F, use option "--ori 2" 

  * **fastq.gz files:** the sequencing data files listed in the sample file.  

2. **<u>Run command:</u>** 

   A typical command for paired end data:

   ```
   amplicon.py -v 20m -j 40 -s sampleFile -k keyFile  -o outputDir -c 2 -n 20 -a 0.01 -m 50 -l 100 -d 0 -z 0.1 -r 10
   ```

   A typical command for single end data:

   ```
   amplicon.py -v 20m -j 40 -s sampleFile -k keyFile -o outputDir  -c 2 -n 20 -a 0.01 -m 50 -l 30 -d 0 -z 0.1 -r 10  --mode 3 --ori 2
   ```

   

   Using "-g" and "-w" parameters to control allele discovery and allele calling mode

   ```
   #"-g" not set.
   If "-g" not set, "-w" will be ignored. With no known allele sequences availabe, the software will call alleles starting from allele ID #1 for each locus.
   
   #"-w 0 -g known_HaplotypeAllele.fasta"
   Call alleles from input file known_HaplotypeAllele.fasta. Alleles that are not present in the fasta file will not be called.
   
   #"-w 1 -g known_HaplotypeAllele.fasta"
   Call alleles from input file known_HaplotypeAllele.fasta. Alleles that are not present in the fasta file will be called, and novel alleles will have allele ID starting from #100001
   
   #"-w 2 -g known_HaplotypeAllele.fasta"
   Call alleles from input file known_HaplotypeAllele.fasta. Alleles that are not present in the fasta file will be called, and novel alleles will have allele ID continue from existing allele ID.
   ```

   The script add_allele.py can be used to add novel alleles from a new HaplotypeAllele.fasta to an existing HaplotypeAllele.fasta file. The allele IDs in the old file will be preserved.  Alleles only present in the new file will be added  in the output merged File, with ID of the novel allele continuing from existing allele ID.  After you have the merged file, you will need to recall genotypes with "-g mergedFile.fasta -w 0" to make sure using the right allele ID (You can use " -i 1" option  when recall on same sample data files. It will skip the most time consuming step 1 of the script).

   ```
   add_allele.py -n newFilePath -o oldFilePath -m mergedFilePath
   ```

   

3. **<u>Output files</u>**

   * hap_genotype: A matrix with all genotypes. Each row is a marker, each column is a sample.

   * hap_genotype_matrix: A file transposed from hap_genotype file, and without read count per allele information. 

   * HaplotypeAllele.fasta: A fasta file with haplotype sequence per allele.

     

4. **<u>Parameters for amplicon.py</u>**

  * -h	show this help message and exit
  * -k         key file
  * -s         sample file
  * -i	Skip steps. e.g. "-i 12" to skip steps 1 and 2. the steps are: 1. split reads by primers; 2. identify haplotypes across population, and optionally run PCR error correction if set "-e 1"; 3 call genotypes. Step 1 is the most time consuming step.  You can adjust most of the following parameters and rerun the script by skipping step 1. 
  * -j	Number of simultaneous jobs. Default:8
  * -t	Number of threads per job. Default:1
  * -c	Minimum number of samples per haplotypes. Default:10
  * -n	Maximum number of unique haplotypes per sample to be kept in the first step. Default:20
  * -m	Maximum number of haplotypes per marker in the population. Default:1000
  * -a	Minimum minor allele frequency Default:0.01
  * -l	Minimum haplotype length (after removing the primers. It must be an integer 1 or larger.) Default:50
  * -d	Whether to merge the duplicate samples (identical names in column 1 of the sample file). 1: Merge; 0: Do not merge and the duplicated sample will be named <sampleName>__<index starting from 1>. Default:1
  * -p	Ploidy, default 2
  * -r	Maximum read count ratio between the two alleles in each sample, default 20
  * -g         Set tag fasta file, so known allele ID will be used. 
  * -w       Novel tag mode. 0: no novel alleles; 1 : novel allele ID start from 100001;  2: novel allelele ID continue from existing allele ID. If "-g" not set, novel allele ID start from 1. 
  * -z	Mismatch rate between pcr primer and reads, default 0.1
  * -u       Set to 1 to restart from crashed point in step 1. 
  * -v        Maximum reads to process per sample.  Default -1 to use all reads in the file.
  * --mode   Set input data mode to 1, 2, or 3. 1:paired end fastq.gz file; 2: genome contigs. 3: single end fastq file. Default:1
  * --ori         Set primer sequence orientation. 1: F-R (default); 2: F-F



### 2. Replace sample names in the  hap_genotype file

```
changeName.py -i PATH_of_input_hap_genotype_file -o PATH_of_outputput_hap_genotype_file -n namefile
```

The name file is a tab delimited text file with two columns (no header line): old name and new name.



### 3. Slice & Merge the hap_genotype file

###### slice.py script

You can slice by either a list of plates or a list of samples. This code will slice out samples from the genotyping matrix. 

To slice by a list of individuals, prepare a sample file, with one sample name per line. 

```
slice.py -i inputDirectory -o outputDirectory -s sampleFileName -a 0.001
```

(Use -a to specify minor allele frequency cutoff. Markers with MAF smaller than cutoff will be discarded. Set "-a 0" to keep all markers)

The inputDirectory should contain at least one file named hap_genotype, and optionally markerToSampleReadCountMatrix from the amplicon.py or merge.py script. The sliced file will be written into the outputDirectory.

In the output hap_genotype file, the samples will be ordered based on the provided sample file. 



###### merge.py script

To merge multiple hap_genotype files, run the command: 

```
merge.py -i inputfile1,inputfile2,inputfile3 -o outputfile 
```

(The input files must be the hap_genotype files from amplicon.py. If there are duplicated sample names, only first ones are kept)



**For Vitisgen project,** 

```
To slice by Vitisgen samples with "plateName_wellName", make a file with "plateName_wellName" for each sample, e.g. vDNAcad794B07_E03. Run command: slice.py -i inputDirectory -o outputDirectory -f vitisgenSampleFileName -a 0.01

To slice by Vitisgen plates, prepare a plate name file, with one plate name per line.
Run command: slice.py -i inputDirectory -o outputDirectory -p plateFileName -a MAF 
```



### 4. Convert to VCF file

##### 4.1. Call SNP and INDELs into a VCF file

see https://bitbucket.org/cornell_bioinformatics/amplicon/src/master/hap2realVCF/

4.2. Make a haplotype VCF file

###### to_lep_map.pl script

Lep-MAP3 requires a vcf file as input. This script converts the hap_genotype from previous step to a "haplotype" VCF file, in which up to 4 haplotype alleles are represented with A C G T codes.  There is a lookup table to show the corresponding haplotype allele sequences for "A" "C" "G" and "T".  

  * Run command

    ```
    to_lep_map.pl -g hap_genotype -f minorAlleleFrequency -b blankSample -m maternalSample -p paternalSample -l marker2pos -n familyName  
    ```

  * -b, -m, -p, -l are optional. "-b" "-m" and "-p" are integer index of the blank, maternal and paternal samples in your sample list. If multiple samples, separate the index with comma. The index should be 1-based, so that the first sample is 1. "-l" is to specify the physical positions of the markers. It should be a tab-delimited table with 3 columns: markerName, chr, pos. 



## Authors

* **Qi Sun**
* **Cheng Zou**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments
The development of the software is supported by Vitisgen2 and Vitisgen3. which are funded by USDA NIFA's Specialty Crop Research Initiative.  