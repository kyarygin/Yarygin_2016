from __future__ import print_function
from time import time
import pandas as pd
import numpy as np
import argparse
import sys
import re
import os

def parse_blast_output(blast_filepath, gene_group_type):
    mult = 3. if gene_group_type == 'prot' else 1.
    blast_output_columns = ['qseqid', 'qlen', 'bgi_id', 'slen',
                            'pident', 'length', 'mismatch', 'gapopen',
                            'qstart', 'qend', 'sstart', 'send',
                            'evalue', 'bitscore']

    blast_result = pd.read_csv(blast_filepath,
        names = blast_output_columns,
        sep='\t', header=None
    )
    blast_result = blast_result[(blast_result.pident > 80) &
                                (blast_result.length * mult / blast_result.slen > 0.8) &
                                (blast_result.length * 1. / blast_result.qlen > 0.8) &
                                (blast_result.evalue < 1e-5)]
    return blast_result

def parse_coverage_files(sample_group_folder, group_name, bgi_names):
    samples_names = [name for name in os.listdir(sample_group_folder) if name.endswith('bp_cov.txt')]
    samples_names = [name[:-11] for name in samples_names]
    for sample_name in samples_names:
        sample_coverage_filepath = os.path.join(sample_group_folder, '%s.bp_cov.txt' % sample_name)
        sample_bgi_coverage = pd.read_csv(sample_coverage_filepath,
            sep='\t', header=None, names=['bp_cov']
        )
        sample_bgi_coverage['sum_bp_cov'] = sum(sample_bgi_coverage['bp_cov'])
        sample_bgi_coverage['sample'] = sample_name
        sample_bgi_coverage['bgi_id'] = bgi_names
        sample_bgi_coverage['group_name'] = group_name
        yield sample_bgi_coverage

def create_coverage_file(blast_result_path, BGI_coverage_path, gene_group_name, input_type):
    bgi_names = pd.read_csv('metadata/BGIGeneSet2010_genes.txt', sep='\t', header=None, names=['bgi_id'])
    sample_groups = os.listdir(BGI_coverage_path)

    blast_result = parse_blast_output(
        blast_filepath=blast_result_path,
        gene_group_type=input_type
    )

    headers = ['qseqid', 'sample', 'group_name', 'abund']
    with open('gene_groups_abund/%s.tsv' % gene_group_name, 'w') as f:
        f.write('\t'.join(headers) + '\n')

    for group_name in sample_groups:
        print(group_name)
        coverage_files = parse_coverage_files(
            sample_group_folder=os.path.join(BGI_coverage_path, group_name),
            group_name=group_name,
            bgi_names=bgi_names
        )

        for sample_bgi_coverage in coverage_files:
            coverage = pd.merge(blast_result, sample_bgi_coverage, on='bgi_id')
            coverage['abund'] = coverage['bp_cov'] / coverage['slen'] / coverage['sum_bp_cov']
            coverage = coverage.loc[:, ['qseqid', 'sample', 'group_name', 'abund']]
            coverage = coverage.groupby(['qseqid', 'sample', 'group_name']).sum()
            coverage.reset_index(inplace=True)
            coverage.to_csv('gene_groups_abund/%s.tsv' % gene_group_name,
                                  sep='\t', mode='a',
                                  header=False, index=False)

def blast(input_file, input_type, gene_group_name, n_threads, config_pathes, blastdb_path):
    if not os.path.isfile(input_file):
        print("Error: file %s doesn't exist" % input_file)
        sys.exit()

    if (input_type == 'prot'):
        blast_type = config_pathes['tblastn_path']
    elif (input_type == 'nucl'):
        blast_type = config_pathes['blastn_path']
    else:
        print('Error: incorrect input type: %s' % input_type)
        sys.exit()

    blast_result_path = 'BLAST_result_%s.txt' % gene_group_name
    cmd = "%s -num_threads %d -db %s -query %s -out %s -evalue 1e-5 -outfmt '6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore'" % (blast_type, n_threads, blastdb_path, input_file, blast_result_path)
    os.system(cmd)

    with open('metadata/gene_groups_types.txt', 'a') as f:
        f.write('%s\t%s\n' % (gene_group_name, input_type))

    return blast_result_path

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', type=str,
                        help='input fasta file')
    parser.add_argument('-n', '--n_threads', type=int,
                        help='number of BLAST threads (default: 1)', default=1)
    parser.add_argument('input_type', type=str,
                        help='type of fasta file (prot/nucl)')
    parser.add_argument('gene_group_name', type=str,
                        help='name of gene group')
    parser.add_argument('output_file', type=str,
                        help='output tsv file')
    parser.add_argument('blastdb_path', type=str,
                        help='path to BGI blastdb')
    parser.add_argument('BGI_coverage_path', type=str,
                        help='path to metagenome coverage')
    args = vars(parser.parse_args())

    with open('config') as config:
        config_pathes = dict([line.strip().split('\t') for line in config.readlines()])

    input_file = args['input_file']
    input_type = args['input_type']
    output_file = args['output_file']
    n_threads = args['n_threads']
    gene_group_name = args['gene_group_name']
    blastdb_path = args['blastdb_path']
    BGI_coverage_path = args['BGI_coverage_path']

    blast_result_path = blast(
        input_file=input_file,
        input_type=input_type,
        gene_group_name=gene_group_name,
        n_threads=n_threads,
        config_pathes=config_pathes,
        blastdb_path=blastdb_path
    )

    create_coverage_file(
        blast_result_path=blast_result_path,
        BGI_coverage_path=BGI_coverage_path,
        gene_group_name=gene_group_name,
        input_type=input_type
    )


