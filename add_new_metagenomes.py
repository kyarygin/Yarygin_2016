import argparse
import sys
import os


if __name__ == '__main__':

    with open('config') as config:
        config_pathes = dict([line.strip().split('\t') for line in config.readlines()])

    samtools_path = config_pathes['samtools_path']
    bowtie_path = config_pathes['bowtie_path']
    genome_coverage_bed_path = config_pathes['genome_coverage_bed_path']
    get_cov_bed_hist_path = config_pathes['get_cov_bed_hist_path']


    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--n_threads', type=int,
                        help='number of bowtie threads (default: 20)', default=20)
    parser.add_argument('group_name', type=str,
                        help='matagenomes group name')
    parser.add_argument('read_file', metavar='read_file', type=str, nargs='+',
                        help='read files')
    args = vars(parser.parse_args())

    # print args

    group_name = args['group_name']
    read_files = args['read_file']
    n_threads =  args['n_threads']

    print group_name
    print read_files
    print n_threads


    os.mkdir('BGI_coverage/%s' % group_name)

    with open('./metadata/BGIGeneSet2010_genes.txt') as f:
        bgi_names = [name.strip() for name in f]


    for reads_path in read_files[:]:

        basename = os.path.basename(reads_path)
        name, ext = os.path.splitext(basename)

        temp_folder = './temp_%s/' % name
        os.mkdir(temp_folder)

        if ext == '.fq' or ext == '.fastq':
            key = '-q'
        elif ext == '.csfasta':
            key = '-f -C'
        else:
            key = '-f'

        if ext == 'csfasta':
            index_path = './index/BGIGeneSet2010_color'
        else:
            index_path = './index/BGIGeneSet2010'

        sam_file = os.path.join(temp_folder, '%s.sam' % name)
        bam_file = os.path.join(temp_folder, '%s.bam' % name)
        log_file = os.path.join(temp_folder, '%s.log' % name)
        coverage_file = os.path.join(temp_folder, '%s.coverage' % name)
        output_file_prefix = os.path.join(temp_folder, '%s.sorted' % name)
        bp_cov_file = os.path.join(temp_folder, '%s.bp_cov.txt' % name)
        pos_cov_file = os.path.join(temp_folder, '%s.pos_cov.txt' % name)

        map_cmd = '%s %s -S -t -v 3 -k 1 --threads %d %s %s %s 2>%s' % \
                         (bowtie_path,
                          key,
                          n_threads,
                          index_path,
                          reads_path,
                          sam_file,
                          log_file)
        os.system(map_cmd)

        sam_to_bam_cmd = '%s import ./index/BGI_GeneSet20090523.fa.fai %s %s' % \
                         (samtools_path,
                          sam_file,
                          bam_file)
        os.system(sam_to_bam_cmd)

        sort_cmd = '%s sort %s %s' % \
                   (samtools_path,
                    bam_file,
                    output_file_prefix)
        os.system(sort_cmd)

        bam_file = output_file_prefix + '.bam'

        coverage_cmd = '%s -ibam %s -g ./index/BGI_GeneSet20090523.fa.genome > %s' % \
                       (genome_coverage_bed_path,
                        bam_file,
                        coverage_file)
        os.system(coverage_cmd)

        bp_pos_cmd = 'python %s %s ./index/BGI_GeneSet20090523.fa.genome %s %s' % \
                     (get_cov_bed_hist_path,
                      coverage_file,
                      bp_cov_file,
                      pos_cov_file)
        os.system(bp_pos_cmd)

        with open(bp_cov_file) as f:
            bgi_abund = dict([x.split('\t') for x in f.readlines()])

        with open('BGI_coverage/%s/%s.txt' % (group_name, name), 'w') as out:
            for bgi_name in bgi_names:
                out.write(bgi_abund[bgi_name])

        os.system('rm -r %s' % temp_folder)
