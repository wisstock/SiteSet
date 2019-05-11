#!/usr/bin/env python3

import sys
import subprocess
import csv
import re
import glob
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord


# Splice sites (ss) extractions, pars by MaxEnt algorithm and CSV dataframe generating
# Copyright (C) 2018 Borys Olifirov
# Institure of Molecular Biology and Genetics of NASU, Systems Biology Research Group
# 
# Variables description:
# f_... or ..._f_... - 5' ss data
# t_... or ..._t_... - 3' ss data
#
# V0.1 - complete splice sites extraction function and correct dataframe generating 


def ss_extract(genbank_in):  # Extraction of list of individual 5' and 3' ss
	seq_record = SeqIO.read(genbank_in, 'genbank')
	f_loc = []
	t_loc = []
	f_list = []
	t_list = []
	f_ss_seq = []
	t_ss_seq = []

	# Loop over the gen file, get the mRNA starts and ends position for '+' strand
	for feature in seq_record.features:
		if feature.type == 'mRNA':
			for exon_location in feature.location.parts:
				t_loc.append(int(exon_location.start))
				f_loc.append(int(exon_location.end))
			break

	with open('seq_log', 'a') as file:  # Adding type, genID and DBref of selected feature to log file
				file.write(genbank_in.split('.')[0])
				file.write(' ' + feature.type + ' ')
				file.write(str(feature.qualifiers.get("db_xref")) + '\n')
				file.write('5` intron end location:' + str(f_loc) + '\n')
				file.write('3` intron end location:' + str(t_loc) + '\n \n')

	# 5' ss extraction
	for n in f_loc:
		if n > 15 and n < len(seq_record.seq):
			f_list.append(SeqFeature.FeatureLocation(n-3, n+6))
	for loc in f_list:
		f_ss_seq.append(str(loc.extract(seq_record.seq)))

	# 3' ss extraction
	for n in t_loc:
		if n > 15:
			t_list.append(SeqFeature.FeatureLocation(n-20, n+3))
	for loc in t_list:
		t_ss_seq.append(str(loc.extract(seq_record.seq)))


	return f_ss_seq, t_ss_seq


def maxent_calc(f_ss, t_ss):  # Collecting results of MaxEnt perl script score calculating
	f_file = open('f_result_seq', 'w')
	f_file.close()
	t_file = open('t_result_seq', 'w')
	t_file.close()
	f_score = []
	t_score = []

	# SU-KA this fucking script needs file but not string as input data
	for one_of_f in f_ss:  # Generating results file of 5'ss
		f_file = open('f_result_seq', 'a')
		f_file.write(one_of_f + '\n')
		f_file.close()

	f_score = subprocess.run(['perl', 'score5.pl', 'f_result_seq'], stdout = subprocess.PIPE).stdout.decode('utf-8')  # Writing data in string format
	f_score = re.findall(r"[-+]?\d*\.\d+|\d+", f_score)  # Extract float values

	for one_of_t in t_ss:  # Generating results file of 3'ss
		t_file = open('t_result_seq', 'a')
		t_file.write(one_of_t + '\n')
		t_file.close()

	t_score = subprocess.run(['perl', 'score3.pl', 't_result_seq'], stdout = subprocess.PIPE).stdout.decode('utf-8')  # Writing data in string format
	t_score = re.findall(r"[-+]?\d*\.\d+|\d+", t_score) # Extract float values


	return(f_score, t_score)


def df_gen(f_data, t_data, gene_title, f_seqs, t_seqs, df_title):  # Dataframe generating
	n = 0

	# Add 5' ss results to csv file
	for one_f_score in f_data:
		one_f_row = [one_f_score, n + 1, '5', f_seqs[n], gene_title.split('.')[0]]
		with open(df_title, 'a', newline = '') as res_file:
			writer = csv.writer(res_file)
			writer.writerow(one_f_row)
		n += 1

	n = 0

	# Add 3' ss results to csv file
	for one_t_score in t_data:
		one_t_row = [one_t_score, n + 1, '3', t_seqs[n], gene_title.split('.')[0]]
		with open(df_title, 'a', newline = '') as res_file:
			writer = csv.writer(res_file)
			writer.writerow(one_t_row)
		n += 1

	


df_name = 'ss_df.csv'

with open('seq_log', 'w') as file:  # Init log file
		file.write('Used features details and secondary info: \n \n')

with open(df_name, 'w', newline = '') as res_file:  # Init total dataframe
	writer = csv.writer(res_file)
	writer.writerow(['score', 'int_num', 'end', 'site_seq', 'gene'])

print('Genes progress:')
for file in glob.glob('*.gb'):
	print(file.split('.')[0])
	test_record = file
	f_ss_sequences, t_ss_sequences = ss_extract(test_record)
	f_ss_scores, t_ss_scores = maxent_calc(f_ss_sequences, t_ss_sequences)
	df_gen(f_ss_scores, t_ss_scores, test_record, f_ss_sequences, t_ss_sequences, df_name)


# That's all!