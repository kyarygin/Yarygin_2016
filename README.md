# Introduction
Pipeline for the quantitative profiling of user defined groups of genes in human gut metagenomes. This method is based on quick analysis of a gene coverage matrix obtained by pre-mapping the metagenomic reads to a global gut microbial catalogue.

# Dependencies
Pipeline relies on **blastn**, **tblastn**, **bowtie**, **samtools** and **genomeCoverageBed** from bedtools inside. If executables are not in your $PATH, please, specify full pathes in `config.json`

# Quick Guide
Download or clone the files from repository and run script `download.py` with Python. This will download all required metadata for pipeline usage.

#### Add new metagenomes
Metadata contains precomputed abundances for each group of metagenomes described in the article. To add new group you need to run `get_genes_abund.py`:
```
$ python add_new_metagenomes.py group_name sample_1.fastq sample_2.fastq sample_3.fastq
```

#### Calculate gene abundances
To calculate the abundances of your gene group run
```
$ python get_genes_abund.py genes.fasta nucl
```
You can specify groups of metagenomes with `-mg` or `--metagen-group` keys. Fow example:
```
$ python get_genes_abund.py genes.fasta nucl -mg USA RUS DEN_control
```
You can see all available metagenome groups by executing: `python get_genes_abund.py -h`. By default all groups of the metagenomes described in article are available.

#### Calculate discrepancy ratios
To calculate the discrepancy ratio for your gene group run
```
$ python get_discrepancy_ratio.py genes.fasta sample_1.fastq sample_2.fastq sample_3.fastq
```
This will produce csv-file with discrepancy ratios for your gene group in the selected samples. Low ratio (<0.1) means that gene catalog is sufficiently representative for this gene group. Discrepancy ratio close to 1 means that the gene group is underrepresented in the catalog and the results should be interpreted with caution.

# Other
* Use key `-bt2` if you want to use Bowtie2 instead of Bowtie
* The `fasta` folder contains the sequences for all gene groups described in the article in FASTA format
* The `gene_groups_abund` folder contains gene abundance tables for all gene groups in all samples used in study. Each file contains abundance of a single gene group in all metagenomes. The table consists of 4 columns: gene name, metagenome name, group of metagenomes name and relative abundance value.
* Run `python get_genes_abund.py -h` and `python add_new_metagenomes.py -h` to get help

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
                        number of Bowtie threads (default: 20)
  -bt2                  use Bowtie2 instead of Bowtie

```
**get_genes_abund.py**
```
$ python get_genes_abund.py [-h] [-o OUTPUT_FOLDER] [-n N_THREADS]
                            [-mg METAGEN_GROUP [METAGEN_GROUP ...]]
                            input_file input_type

positional arguments:
  input_file            path to input FASTA file
  input_type            type of FASTA file (prot/nucl)

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
  gene_group_file       path to input FASTA file
  read_file             path to read files

optional arguments:
  -h, --help            show this help message and exit
  -n N_THREADS, --n_threads N_THREADS
                        number of Bowtie threads (default: 20)
  -o OUTPUT_FOLDER, --output-folder OUTPUT_FOLDER
                        path to output folder (default: current dir)
  -bt2                  use Bowtie2 instead of Bowtie

```
