### hap2realVCF

Convert the amplicon.py output to real VCF. 

##### Dependencies:

- PYTHON modules: screed, scikit-allele, pandas, scipy
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

Create vcf file, input files including: the directory from hap2reads.py, and bwa indexed reference genome

```
hap2vcf.sh path_of_readOutputDir full_path_of_mygenome.fasta number_of_parallel_jobs
```

The output vcf file is unfiltered.vcf. 

If the vcf file is ok, you want to the delete the fastq and directories to save storage space. The bam directory should be named "bam" in current directory

##### Step 3: filter the vcf file

a. filter sites by missing data. (in the example, keep sites with missing data<0.5)

```
bcftools filter -i 'F_MISSING<0.5'  -O v -o filtered_sites.vcf unfiltered.vcf
```

b. filter samples by missing data (in the example, keep samples with missing data<0.5)

```
vcftools --missing-indv --vcf filtered_sites.vcf

awk '{if ($5<0.50) print}' out.imiss |cut -f1 > keep

bcftools view -S keep filtered_sites.vcf > filtered_sites_sample.vcf
```

The filtered vcf file is filtered_sites_sample.vcf 

c. check minor allele frequency per site (optional)

```
bcftools +fill-tags filtered_sites_sample.vcf > final.vcf

bcftools query -f '%CHROM\t%POS\t%AF\t%NS\n'  final.vcf > maf.stat.txt
```

In the maf.stat.txt, the first column is the allele frequency, NS is number of samples with data 
