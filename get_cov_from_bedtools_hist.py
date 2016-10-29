# -*- coding: utf-8 -*-
#!/usr/bin/python

#################################
# Get total coverage (bp_coverage) and number of covered nucleotides (pos_coverage) from bedtools coverage histogram.
# The resulting coverage is sorted alphabetically by seq name.
# Zero-covered sequences are included.
#################################
import sys
import string

file_cov_hist = sys.argv[1]
file_bedtools_genome = sys.argv[2]
file_out_bp_cov = sys.argv[3]
file_out_pos_cov = sys.argv[4]

list_bp_cov=dict()
list_pos_cov=dict()

# get all sequences names
h_file_bedtools_genome = open(file_bedtools_genome, 'rU')
for pline in h_file_bedtools_genome:
	words=pline.split()
	list_bp_cov[words[0]] = 0
	list_pos_cov[words[0]] = 0

# get coverages for seqs
h_file_cov_hist = open(file_cov_hist, 'rU')
for pline in h_file_cov_hist:
	words=pline.split()
	if words[0] != 'genome': # special reserved word of bedtools (for summary coverage)
		list_bp_cov[words[0]] += int(words[1])*int(words[2])		
		if words[1] != '0':
			list_pos_cov[words[0]] += int(words[2])

h_file_out_bp_cov = open(file_out_bp_cov, 'w')
# zero-covered seqs are omitted in histogram!
# so fill the lacking and output
for key in sorted(list_bp_cov.keys()):
	h_file_out_bp_cov.write(key + '\t' + str(list_bp_cov[key]) + '\n')
h_file_out_bp_cov.close()

# get pos coverage
h_file_out_pos_cov = open(file_out_pos_cov, 'w')
for key in sorted(list_pos_cov.keys()):
	h_file_out_pos_cov.write(key + '\t' + str(list_pos_cov[key]) + '\n')
h_file_out_pos_cov.close()

