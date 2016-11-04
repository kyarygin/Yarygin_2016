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

def create_coverage_file(blast_result_path, bgi_coverage_path, gene_group_name, output_folder, input_type, mg_groups):
    bgi_names_path = os.path.join('metadata', 'BGIGeneSet2010_genes.txt')
    bgi_names = pd.read_csv(bgi_names_path, sep='\t', header=None, names=['bgi_id'])

    blast_result = parse_blast_output(
        blast_filepath=blast_result_path,
        gene_group_type=input_type
    )

    headers = ['qseqid', 'sample', 'group_name', 'abund']
    gene_group_abund_path = os.path.join(output_folder, '{}.tsv'.format(gene_group_name))
    with open(gene_group_abund_path, 'w') as f:
        f.write('\t'.join(headers) + '\n')

    for group_name in mg_groups:
        coverage_files = parse_coverage_files(
            sample_group_folder=os.path.join(bgi_coverage_path, group_name),
            group_name=group_name,
            bgi_names=bgi_names
        )

        for sample_bgi_coverage in coverage_files:
            coverage = pd.merge(blast_result, sample_bgi_coverage, on='bgi_id')
            coverage['abund'] = coverage['bp_cov'] / coverage['slen'] / coverage['sum_bp_cov']
            coverage = coverage.loc[:, ['qseqid', 'sample', 'group_name', 'abund']]
            coverage = coverage.groupby(['qseqid', 'sample', 'group_name']).sum()
            coverage.reset_index(inplace=True)
            coverage.to_csv(gene_group_abund_path,
                            sep='\t', mode='a',
                            header=False, index=False)

def blast(input_file, input_type, gene_group_name, n_threads, config_pathes, blastdb_path):
    if not os.path.isfile(input_file):
        sys.stdout.write('Error: file {} doesn\'t exist\n'.format(input_file))
        sys.exit()

    if (input_type == 'prot'):
        blast_type = config_pathes['tblastn_path']
    elif (input_type == 'nucl'):
        blast_type = config_pathes['blastn_path']
    else:
        sys.stdout.write('Error: incorrect input type: {}\n'.format(input_type))
        sys.exit()

    blast_result_path = os.path.join('blast_results', 'BLAST_result_{}.txt'.format(gene_group_name))
    outfmt = '\'6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore\''
    cmd = '{blast_type} -num_threads {n_threads} -db {blastdb_path} -query {input_file} -out {blast_result_path} -evalue 1e-5 -outfmt {outfmt}'.format(
        blast_type=blast_type,
        n_threads=n_threads,
        blastdb_path=blastdb_path,
        input_file=input_file,
        blast_result_path=blast_result_path,
        outfmt=outfmt
    )
    os.system(cmd)

    return blast_result_path

def parse_arguments():
    available_groups = set(os.listdir('BGI_coverage'))
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input_file', type=str,
                        help='path to input fasta file')
    parser.add_argument('input_type', type=str,
                        help='type of fasta file (prot/nucl)')
    parser.add_argument('-o', '--output-folder', type=str, default='.',
                        help='path to output folder (default: current dir)')
    parser.add_argument('-n', '--n-threads', type=int,
                        help='number of BLAST threads (default: 1)', default=1)
    mg_help = ['select groups of metagenomes',
               'available: {}'.format(', '.join('\'{}\''.format(group_name) for group_name in available_groups)),
               'default: all']
    parser.add_argument('-mg', '--metagen-group', nargs='+',
                        default=available_groups,
                        help='\n'.join(mg_help))
    args = vars(parser.parse_args())

    if not os.path.exists(args['input_file']):
        raise Exception('Input file doesn\'t exist')

    if not os.path.exists(args['output_folder']) or not os.path.isdir(args['output_folder']):
        raise Exception('Output file folder doesn\'t exist')

    if args['input_type'] not in {'prot', 'nucl'}:
        raise Exception('Invalid input type: {}'.format(args['input_type']))

    if args['n_threads'] <= 0:
        raise Exception('Invalid number of threads: {}'.format(args['n_threads']))

    wrong_groups = [group_name for group_name in args['metagen_group'] if group_name not in available_groups]
    if wrong_groups:
        raise Exception('Invalid groups of metagenomes: {}'.format(', '.join(wrong_groups)))

    return args

if __name__ == '__main__':
    args = parse_arguments()
    mg_groups = args['metagen_group']
    input_type = args['input_type']
    input_file = args['input_file']
    output_folder = args['output_folder']
    n_threads = args['n_threads']

    gene_group_name = os.path.splitext(os.path.basename(input_file))[0]

    blastdb_path = os.path.join('BGI_blastdb', 'BGIGeneSet2010')
    bgi_coverage_path = 'BGI_coverage'

    with open('config.json') as f:
        config_pathes = eval(''.join(f.readlines()))

    sys.stdout.write('Running BLAST ... \n')
    blast_result_path = blast(
        input_file=input_file,
        input_type=input_type,
        gene_group_name=gene_group_name,
        n_threads=n_threads,
        config_pathes=config_pathes,
        blastdb_path=blastdb_path
    )

    sys.stdout.write('Creating coverage file ... \n')
    create_coverage_file(
        blast_result_path=blast_result_path,
        bgi_coverage_path=bgi_coverage_path,
        gene_group_name=gene_group_name,
        output_folder=output_folder,
        input_type=input_type,
        mg_groups=mg_groups
    )
