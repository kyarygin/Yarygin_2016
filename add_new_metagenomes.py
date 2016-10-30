from __future__ import print_function
import argparse
import sys
import os

def delete_folder(top):
    for root, dirs, files in os.walk(top, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))
    os.rmdir(top)

if __name__ == '__main__':

    with open('config.json') as f:
        config_pathes = eval(''.join(f.readlines()))

    samtools_path = config_pathes['samtools_path']
    bowtie_path = config_pathes['bowtie_path']
    genome_coverage_bed_path = config_pathes['genome_coverage_bed_path']
    get_cov_bed_hist_path = 'get_cov_from_bedtools_hist.py'


    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--n_threads', type=int,
                        help='number of bowtie threads (default: 20)', default=20)
    parser.add_argument('group_name', type=str,
                        help='matagenomes group name')
    parser.add_argument('read_file', metavar='read_file', type=str, nargs='+',
                        help='read files')
    args = vars(parser.parse_args())

    group_name = args['group_name']
    read_files = args['read_file']
    n_threads =  args['n_threads']

    os.mkdir(os.path.join('BGI_coverage', group_name))

    bgi_names_path = os.path.join('metadata', 'BGIGeneSet2010_genes.txt')
    with open(bgi_names_path) as f:
        bgi_names = [name.strip() for name in f]


    for reads_path in read_files:
        basename = os.path.basename(reads_path)
        name, ext = os.path.splitext(basename)

        sys.stdout.write('Processing {} ... \n'.format(basename))

        temp_folder = 'temp_{}'.format(name)
        os.mkdir(temp_folder)

        if ext == '.fq' or ext == '.fastq':
            key = '-q'
        elif ext == '.csfasta':
            key = '-f -C'
        else:
            key = '-f'

        if ext == 'csfasta':
            index_path = os.path.join('index', 'BGIGeneSet2010_color')
        else:
            index_path = os.path.join('index', 'BGIGeneSet2010')

        sam_file = os.path.join(temp_folder, '{}.sam'.format(name))
        bam_file = os.path.join(temp_folder, '{}.bam'.format(name))
        log_file = os.path.join(temp_folder, '{}.log'.format(name))
        coverage_file = os.path.join(temp_folder, '{}.coverage'.format(name))
        output_file_prefix = os.path.join(temp_folder, '{}.sorted'.format(name))
        bp_cov_file = os.path.join(temp_folder, '{}.bp_cov.txt'.format(name))
        pos_cov_file = os.path.join(temp_folder, '{}.pos_cov.txt'.format(name))
        bgi_fai_file = os.path.join('index', 'BGIGeneSet2010.fasta.fai')
        bgi_genome_file = os.path.join('index', 'BGIGeneSet2010.fasta.genome')

        sys.stdout.write('Mapping on catalog ... \n')
        map_cmd = '{bowtie_path} {key} -S -t -v 3 -k 1 --threads {n_threads} {index_path} {reads_path} {sam_file} 2>{log_file}'.format(
            bowtie_path=bowtie_path,
            key=key,
            n_threads=n_threads,
            index_path=index_path,
            reads_path=reads_path,
            sam_file=sam_file,
            log_file=log_file
        )
        os.system(map_cmd)

        sys.stdout.write('Processing .sam file ... \n')
        sam_to_bam_cmd = '{samtools_path} import {bgi_fai_file} {sam_file} {bam_file}'.format(
            samtools_path=samtools_path,
            bgi_fai_file=bgi_fai_file,
            sam_file=sam_file,
            bam_file=bam_file
        )
        os.system(sam_to_bam_cmd)

        sort_cmd = '{samtools_path} sort {bam_file} {output_file_prefix}'.format(
            samtools_path=samtools_path,
            bam_file=bam_file,
            output_file_prefix=output_file_prefix
        )
        os.system(sort_cmd)

        bam_file = output_file_prefix + '.bam'

        sys.stdout.write('Calculating coverage ... \n')
        coverage_cmd = '{genome_coverage_bed_path} -ibam {bam_file} -g {bgi_genome_file} > {coverage_file}'.format(
            genome_coverage_bed_path=genome_coverage_bed_path,
            bam_file=bam_file,
            bgi_genome_file=bgi_genome_file,
            coverage_file=coverage_file
        )
        os.system(coverage_cmd)

        bp_pos_cmd = 'python {get_cov_bed_hist_path} {coverage_file} {bgi_genome_file} {bp_cov_file} {pos_cov_file}'.format(
            get_cov_bed_hist_path=get_cov_bed_hist_path,
            coverage_file=coverage_file,
            bgi_genome_file=bgi_genome_file,
            bp_cov_file=bp_cov_file,
            pos_cov_file=pos_cov_file
        )
        os.system(bp_pos_cmd)

        with open(bp_cov_file) as f:
            bgi_abund = dict([x.split('\t') for x in f.readlines()])

        abund_file_path = os.path.join('BGI_coverage', group_name, '{}.txt'.format(name))
        with open(abund_file_path, 'w') as out:
            for bgi_name in bgi_names:
                out.write(bgi_abund[bgi_name])

        delete_folder(temp_folder)
