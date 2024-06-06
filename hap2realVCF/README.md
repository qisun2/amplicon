### hap2realVCF

Convert the amplicon.py output to real VCF. 

##### Dependencies:

- PYTHON modules: pandas, Biopython
- bwa
- bcftools
- vcftools (optional, for sample missing data statistics)
- gnu parallel

##### Installation
Use pip to install the required python modules. All you need is this python script file.

##### Step 1: create fastq files for each sample, using two output files from amplicon.py: hap_genotype HaplotypeAllele.fasta. The output fastq files will be created in the directory readOutputDir

```
hap2reads.py hap_genotype HaplotypeAllele.fasta readOutputDir
```

##### Step 2: create vcf file from the fastq files

Prepare reference genome (bwa indexed genome) :

```
bwa index mygenome.fasta
```

Create vcf file, input files including: the directory from hap2reads.py, and bwa indexed reference genome. The last parameter is the number of jobs running in parallel. This number is dependent on how many CPU cores on your computer. If not set, default to 4. Each job will use 2 cpu cores. Do not set this number too large. Otherwise you will run into memory problems. (use hap2vcf_mac.sh if you run on a Mac) 

```
hap2vcf.sh path_of_readOutputDir full_path_of_mygenome.fasta number_of_parallel_jobs
```

The output vcf file is unfiltered.vcf. 

If the vcf file looks ok, you might want to delete the fastq and bam directories to save storage space. The fastq directory is the output from hap2read.py, and the bam directory is named "bam" and located in current directory.

##### Step 3: filter the vcf file

a. filter sites by missing data. (in the example, keep sites with missing data<0.3, and maf > 0.01)

```
bcftools view -q 0.01:minor -i 'F_MISSING<0.3' -o filtered.vcf.gz -O z unfiltered.vcf.gz
```

b. filter samples by missing data (in the example, keep samples with missing data<0.5)

```
#get sample list with missing data
bcftools query -f '[%SAMPLE\t%GT\n]' filtered.vcf.gz | \
awk '{count[$1]++; if($2 == "./.") missing[$1]++} END {for (sample in count) print sample, missing[sample]/count[sample]}' > sample_missing_data.txt

#filter
awk '$2 > 0.5 {print $1}' sample_missing_data.txt > samples_to_remove.txt

#filter these samples
bcftools view --samples-file ^samples_to_remove.txt -Oz -o sample_filtered.vcf.gz filtered.vcf.gz

```

The filtered vcf file is filtered_sites_sample.vcf 

c. get stats

```
bcftools +fill-tags -O z -o final.vcf.gz  filtered.vcf.gz

bcftools query -f '%CHROM\t%POS\t%AF\t%NS\n'  final.vcf.gz > maf_ns.stat.txt
```

In the maf_ns.stat.txt, the 3rd column is the allele frequency, 4th column is the number of samples with genotyping data.
