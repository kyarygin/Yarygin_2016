# Introduction
Pipeline for the quantitative profiling of user defined groups of genes in human gut metagenomes. This method is based on the quick analysis of a gene coverage matrix obtained by pre-mapping the metagenomic reads to a global gut microbial catalogue.

# Dependencies
Pipeline relies on **blastn**, **tblastn**, **bowtie**, **samtools** and **genomeCoverageBed** from bedtools inside. If executables not in your $PATH, please, specify full pathes in `config.json`

# Guide
Download or clone files from repository and run script `download.py` with python. This will download all required metadata for pipeline usage.

Metadata contains precomputed abundances for metagenome groups stated in article. To add new group you need to run `get_genes_abund.py`:
```
$ python add_new_metagenomes.py group_name sample_1.fastq sample_2.fastq sample_3.fastq
```

To get abundances of your gene group run
```
$ python get_genes_abund.py genes.fasta nucl
```
You can specify metagenome groups with `-mg` or `--metagenomic-group` keys. Fow example:
```
$ python get_genes_abund.py genes.fasta nucl -mg USA RUS DEN_control
```
You cas see all available metagenome groups executing `python get_genes_abund.py -h`. By default all metagenome groups stated in article are available.
