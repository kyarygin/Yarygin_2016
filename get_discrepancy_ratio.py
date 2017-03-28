import pysam
import argparse
import sys
import re
import os

def delete_folder(top):
    for root, dirs, files in os.walk(top, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))
    os.rmdir(top)

def build_index(input_fasta_path, bowtie_build_path):
    gene_group_name, _ = os.path.splitext(os.path.basename(input_fasta_path))
    index_path = os.path.join('temp', 'index', gene_group_name)
    build_cmd = '{bowtie_build_path} {input_fasta_path} {index_path}'.format(
        bowtie_build_path=bowtie_build_path,
        input_fasta_path=input_fasta_path,
        index_path=index_path
    )
    os.system(build_cmd)
    return index_path

def map_to_bgi(reads_path, bowtie_path, n_threads, bt2):
    name, ext = os.path.splitext(os.path.basename(reads_path))

    if ext == '.fq' or ext == '.fastq':
        key = '-q'
    elif ext == '.csfasta':
        key = '-f -C'
    else:
        key = '-f'

    if bt2:
        index_path = os.path.join('index', 'BGIGeneSet2010_bt2')
    elif ext == 'csfasta':
        index_path = os.path.join('index', 'BGIGeneSet2010_color')
    else:
        index_path = os.path.join('index', 'BGIGeneSet2010')

    sam_file = os.path.join('temp', 'bgi_sam', '{}.sam'.format(name))
    log_file = os.path.join('temp', 'logs', '{}.bgi_mapping.log'.format(name))
    unmapped_fasta = os.path.join('temp', 'unmapped', '{}.unmapped{}'.format(name, ext))

    if bt2:
        map_cmd = '{bowtie_path} -x {index_path} -U {reads_path} -S {sam_file} -k 1 -p {n_threads} --un {unmapped_fasta} 2>{log_file}'.format(
            bowtie_path=bowtie_path,
            index_path=index_path,
            reads_path=reads_path,
            sam_file=sam_file,
            n_threads=n_threads,
            unmapped_fasta=unmapped_fasta,
            log_file=log_file
        )
    else:
        map_cmd = '{bowtie_path} {key} -S -t -v 3 -k 1 --threads {n_threads} --un {unmapped_fasta} {index_path} {reads_path} {sam_file} 2>{log_file}'.format(
            bowtie_path=bowtie_path,
            key=key,
            n_threads=n_threads,
            unmapped_fasta=unmapped_fasta,
            index_path=index_path,
            reads_path=reads_path,
            sam_file=sam_file,
            log_file=log_file
        )
    os.system(map_cmd)

def map_to_raw(reads_path, index_path, bowtie_path, n_threads):
    name, ext = os.path.splitext(os.path.basename(reads_path))

    sam_file = os.path.join('temp', 'raw_sam', '{}.sam'.format(name))
    log_file = os.path.join('temp', 'logs', '{}.raw_mapping.log'.format(name))
    unmapped_fasta = os.path.join('temp', 'unmapped', '{}.unmapped{}'.format(name, ext))

    if ext == '.fq' or ext == '.fastq':
        key = '-q'
    elif ext == '.csfasta':
        key = '-f -C'
    else:
        key = '-f'

    if bt2:
        map_cmd = '{bowtie_path} -x {index_path} -U {reads_path} -S {sam_file} -k 1 -p {n_threads} --no-unal 2>{log_file}'.format(
            bowtie_path=bowtie_path,
            index_path=index_path,
            reads_path=unmapped_fasta,
            sam_file=sam_file,
            n_threads=n_threads,
            log_file=log_file
        )
    else:
        map_cmd = '{bowtie_path} {key} -S -t -v 3 -k 1 --threads {n_threads} {index_path} {reads_path} {sam_file} 2>{log_file}'.format(
            bowtie_path=bowtie_path,
            key=key,
            n_threads=n_threads,
            index_path=index_path,
            reads_path=unmapped_fasta,
            sam_file=sam_file,
            log_file=log_file
        )
    os.system(map_cmd)

def calc_discrepancy_ratio(reads_path):
    name, ext = os.path.splitext(os.path.basename(reads_path))

    bgi_sam_file = os.path.join('temp', 'bgi_sam', '{}.sam'.format(name))
    raw_sam_file = os.path.join('temp', 'raw_sam', '{}.sam'.format(name))

    bgi_read_count = 0
    raw_read_count = 0
    for sam_record in pysam.AlignmentFile(bgi_sam_file):
        if sam_record.reference_id != -1:
            bgi_read_count += 1
    for sam_record in pysam.AlignmentFile(raw_sam_file):
        if sam_record.reference_id != -1:
            raw_read_count += 1
    if bgi_read_count == 0:
        return 'Inf'
    else:
        return 1. * raw_read_count / bgi_read_count


def parse_arguments():
    available_groups = set(os.listdir('BGI_coverage'))
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('gene_group_file', type=str,
                        help='path to input fasta file')
    parser.add_argument('read_file', metavar='read_file', type=str, nargs='+',
                        help='path to read files')
    parser.add_argument('-n', '--n_threads', type=int,
                        help='number of bowtie threads (default: 20)', default=20)
    parser.add_argument('-o', '--output-folder', type=str, default='.',
                        help='path to output folder (default: current dir)')
    parser.add_argument('-bt2', action='store_true',
                        help='use bowtie2 instead of bowtie')

    args = vars(parser.parse_args())
    bt2 = args['bt2']

    if not os.path.exists(args['gene_group_file']):
        raise Exception('Input file doesn\'t exist')

    if not os.path.exists(args['output_folder']) or not os.path.isdir(args['output_folder']):
        raise Exception('Output file folder doesn\'t exist')

    if args['n_threads'] <= 0:
        raise Exception('Invalid number of threads: {}'.format(args['n_threads']))

    return args

if __name__ == '__main__':
    args = parse_arguments()
    gene_group_file = args['gene_group_file']
    read_files = args['read_file']
    output_folder = args['output_folder']
    n_threads = args['n_threads']
    with open('config.json') as f:
        config_pathes = eval(f.read())

    blastdb_path = os.path.join('BGI_blastdb', 'BGIGeneSet2010')
    bgi_coverage_path = 'BGI_coverage'

    os.mkdir('temp')
    os.mkdir(os.path.join('temp', 'unmapped'))
    os.mkdir(os.path.join('temp', 'bgi_sam'))
    os.mkdir(os.path.join('temp', 'raw_sam'))
    os.mkdir(os.path.join('temp', 'logs'))
    os.mkdir(os.path.join('temp', 'index'))

    index_path = build_index(gene_group_file, config_pathes['bowtie-build_path'])
    output_file = os.path.join(output_folder, 'discrepancy_ratio.csv')
    with open(output_file, 'w') as f:
        f.write('readset,discrepancy_ratio\n')

    for read_file in read_files:
        map_to_bgi(reads_path=read_file,
                   bowtie_path=config_pathes['bowtie_path'],
                   n_threads=n_threads,
                   bt2=bt2)
        map_to_raw(reads_path=read_file,
                   index_path=index_path,
                   bowtie_path=config_pathes['bowtie_path'],
                   n_threads=n_threads,
                   bt2=bt2)
        discrepancy_ratio = calc_discrepancy_ratio(read_file)
        with open(output_file, 'a') as f:
            f.write('{},{}\n'.format(read_file, discrepancy_ratio))

    delete_folder('temp')
