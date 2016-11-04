# Introduction
Pipeline for the quantitative profiling of user defined groups of genes in human gut metagenomes. This method is based on the quick analysis of a gene coverage matrix obtained by pre-mapping the metagenomic reads to a global gut microbial catalogue.

# Dependencies
Pipeline relies on **blastn**, **tblastn**, **bowtie**, **samtools** and **genomeCoverageBed** from bedtools inside. If executables not in your $PATH, please, specify full pathes in `config.json`

# Quick Guide
Download or clone files from repository and run script `download.py` with python. This will download all required metadata for pipeline usage.

#### Add new metagenomes
Metadata contains precomputed abundances for each group of metagenomes stated in article. To add new group you need to run `get_genes_abund.py`:
```
$ python add_new_metagenomes.py group_name sample_1.fastq sample_2.fastq sample_3.fastq
```

#### Calculate gene abundances
To get abundances of your gene group run
```
$ python get_genes_abund.py genes.fasta nucl
```
You can specify groups of metagenomes with `-mg` or `--metagenomic-group` keys. Fow example:
```
$ python get_genes_abund.py genes.fasta nucl -mg USA RUS DEN_control
```
You can see all available metagenome groups executing `python get_genes_abund.py -h`. By default all metagenome groups stated in article are available.

# Script usage

```
$ python add_new_metagenomes.py [-h] [-n N_THREADS]
                                group_name read_file [read_file ...]

positional arguments:
  group_name            matagenomes group name
  read_file             read files

optional arguments:
  -h, --help            show help message and exit
  -n N_THREADS, --n_threads N_THREADS
                        number of bowtie threads (default: 20)

```
```
$ python get_genes_abund.py [-h] [-o OUTPUT_FOLDER] [-n N_THREADS]
                            [-mg METAGENOMIC_GROUP [METAGENOMIC_GROUP ...]]
                            input_file input_type

positional arguments:
  input_file            path to input fasta file
  input_type            type of fasta file (prot/nucl)

optional arguments:
  -h, --help            show help message and exit
  -o OUTPUT_FOLDER, --output-folder OUTPUT_FOLDER
                        path to output folder (default: current dir)
  -n N_THREADS, --n-threads N_THREADS
                        number of BLAST threads (default: 1)
  -mg METAGENOMIC_GROUP [METAGENOMIC_GROUP ...], --metagenomic-group METAGENOMIC_GROUP [METAGENOMIC_GROUP ...]
                        select metagenomic groups
                        available: 'USA', 'CHI', 'SPN_CD', 'SPN_control', 'DEN_obese', 'DEN_control', 'RUS_plus', 'RUS', 'SPN_UC'
                        default: all

```

# Other
* `fasta` folder contains all genes used in study in fasta-format
* `gene_groups_abund` folder contains gene abundance tables for all gene groups in all samples used in study
* Run `python get_genes_abund.py -h` and `python add_new_metagenomes.py -h` to get usage help
