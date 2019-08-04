# Analyze Amplicon

The python script (amplicon.py)  analyzes multi-plexed amplicon sequencing data, optimized for high throughput IDT rhAmpSeq data. Another script (to_lep_map.pl) converta the output to a VCF file that can be loaded into LepMap3 software for genetic linkage map construction. The alleles in the VCF file (A, C, G, T) are symbols represent up to 4 haplotytpe alleles. The sequences of the haplotype allleles represented by "ACGT" are provided in the output.

## Getting Started


### Prerequisites
1. Python 3 (tested on python 3.6) and PERL (tested on 5.22)
2. Python modules: 
* BioPython https://biopython.org/
3. Other software. The following commands should in the PATH
* bbmerge.sh: https://sourceforge.net/projects/bbmap/
* cutadapt: https://cutadapt.readthedocs.io/en/stable/
* muscle (optional): https://www.drive5.com/muscle/

### Installing

Install all pre-requisites. Download the amplicon.py and to_lep_map.pl scripts.

### Usage

## Authors

* **Qi Sun**
* **Cheng Zou**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments
