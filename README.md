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
You can specify groups of metagenomes with `-mg` or `--metagen-group` keys. Fow example:
```
$ python get_genes_abund.py genes.fasta nucl -mg USA RUS DEN_control
```
You can see all available metagenome groups executing `python get_genes_abund.py -h`. By default all groups of metagenomes stated in article are available.

#### Calculate discrepancy ratios
To get discrepancy ratio of your gene group run
```
$ python get_discrepancy_ratio.py genes.fasta sample_1.fastq sample_2.fastq sample_3.fastq
```
This will produce csv-file with discrepancy ratios of your gene group in selected samples. Low ratio (<0.1) means that gene catalog is representative to gene group. Discrepancy ratio ~1 means gene group is underrepresented in the catalog.

# Other
* Use key `-bt2` if you want to use bowtie2 instead of bowtie
* `fasta` folder contains all genes used in study in fasta-format
* `gene_groups_abund` folder contains gene abundance tables for all gene groups in all samples used in study. Each file fontains abundace of one gene group in all metagenomes. Table consist of four columns - gene name, metagenome name, group of metagenomes name and relative abundance value.
* Run `python get_genes_abund.py -h` and `python add_new_metagenomes.py -h` to get usage help

# Script usage
**add_new_metagenomes.py**
```
$ python add_new_metagenomes.py [-h] [-n N_THREADS] [-bt2]
                                group_name read_file [read_file ...]

positional arguments:
  group_name            name of new group of metagenomes
  read_file             path to read files

optional arguments:
  -h, --help            show this help message and exit
  -n N_THREADS, --n_threads N_THREADS
                        number of bowtie threads (default: 20)
  -bt2                  use bowtie2 instead of bowtie

```
**get_genes_abund.py**
```
$ python get_genes_abund.py [-h] [-o OUTPUT_FOLDER] [-n N_THREADS]
                            [-mg METAGEN_GROUP [METAGEN_GROUP ...]]
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
  -mg METAGEN_GROUP [METAGEN_GROUP ...], --metagen-group METAGEN_GROUP [METAGEN_GROUP ...]
                        select groups of metagenomes
                        available: 'USA', 'CHI', 'SPN_CD', 'SPN_control', 'DEN_obese', 'DEN_control', 'RUS_plus', 'RUS', 'SPN_UC'
                        default: all

```
**get_discrepancy_ratio.py**
```
python get_discrepancy_ratio.py [-h] [-n N_THREADS] [-o OUTPUT_FOLDER] [-bt2]
                                gene_group_file read_file [read_file ...]

positional arguments:
  gene_group_file       path to input fasta file
  read_file             path to read files

optional arguments:
  -h, --help            show this help message and exit
  -n N_THREADS, --n_threads N_THREADS
                        number of bowtie threads (default: 20)
  -o OUTPUT_FOLDER, --output-folder OUTPUT_FOLDER
                        path to output folder (default: current dir)
  -bt2                  use bowtie2 instead of bowtie

```
