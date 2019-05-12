#!/usr/bin/env python3

import sys
import subprocess
import csv
import re
import glob
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord


"""
Copyright © 2018 Institure of Molecular Biology and Genetics of NASU, Systems Biology Research Group
Copyright © 2018 Borys Olifirov
Author e-mail: omnia.fatum@gmail.com

Splice sites (ss) extractions orient by exon features location,
score calculating by MaxEnt algorithm and output CSV dataframe generating

The script require Python 3.5 or higer, Biopyton 1.68 or higer
and MaxEnt scripts (Yeo G and Burge C.B., Maximum Entropy Modeling
of Short Sequence Motifs with Applications to RNA Splicing Signals,
Journal of Computational Biology, in press; 
Copyright 2004, Gene Yeo, Chris Burge)

Variables description:
f_... or ..._f_...		5' ss data
t_... or ..._t_...		3' ss data

V0.2 - complete splice sites extraction function and correct
       dataframe generating
V0.3   + automatic identification of intron with PAS (in progress)

"""


def ss_extract(genbank_in):
	""" 
	Extraction of list of individual 5' and 3' ss.
	Function takes an input GenBank file name (genbank_in);
	Function returns 5' and 3' ss sequences (f_ss_seq, t_ss_seq)

	"""
	seq_record = SeqIO.read(genbank_in, 'genbank')
	f_loc = []
	t_loc = []
	f_list = []
	t_list = []
	f_ss_seq = []
	t_ss_seq = []
	
	# Loop over the gen file, get the mRNA starts
	# and ends position for '+' strand
	for feature in seq_record.features:
		if feature.type == 'mRNA':
			for exon_location in feature.location.parts:
				t_loc.append(int(exon_location.start))
				f_loc.append(int(exon_location.end))
			break

	with open('seq_log', 'a') as file:  # Adding data to log file
				file.write(genbank_in.split('.')[0])  							# Gene name
				file.write(' ' + feature.type + ' ') 							# Selected feature type
				file.write(str(feature.qualifiers.get("db_xref")) + '\n') 		# GeneID and DBref feature
				file.write('5` intron end location:' + str(f_loc) + '\n')		# 5' intron ends location
				file.write('3` intron end location:' + str(t_loc) + '\n \n')	# 3' intron ends location

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


def maxent_calc(f_ss, t_ss, files_title):  
	""" 
	Collecting results of MaxEnt perl script score calculating.
	Function takes an input lists of 5' and 3' ss sequences (f_ss, t_ss)
	and input GenBank file name (gene_title);
	Function return results of MaxEnt perl script
	in a list format (f_score, t_score) 

	"""
	f_file_title = str(files_title.split('.')[0] + '_f_seq')
	t_file_title = str(files_title.split('.')[0] + '_t_seq')

	with open(f_file_title, 'w') as f_file:
		f_file.write('> ' + files_title.split('.')[0] + '5` ss sequences \n')
	with open(t_file_title, 'w') as t_file:
		t_file.write('> ' + files_title.split('.')[0] + '3` ss sequences \n')

	f_score = []
	t_score = []

	# SU-KA this fucking script needs file but not string as input data
	# with open(f_file_title, 'a') as f_file:  # Generating results file of 5'ss
	# 	f_file.write(one_of_f + '\n')
	for one_of_f in f_ss: 
		f_file = open(f_file_title, 'a')
		f_file.write(one_of_f + '\n')
		f_file.close()

	f_score = subprocess.run(['perl', 'score5.pl', f_file_title], stdout = subprocess.PIPE).stdout.decode('utf-8')  # Writing data in string format
	f_score = re.findall(r"[-+]?\d*\.\d+|\d+", f_score)  # Extract float values

	for one_of_t in t_ss:  # Generating results file of 3'ss
		t_file = open(t_file_title, 'a')
		t_file.write(one_of_t + '\n')
		t_file.close()

	t_score = subprocess.run(['perl', 'score3.pl', t_file_title], stdout = subprocess.PIPE).stdout.decode('utf-8')  # Writing data in string format
	t_score = re.findall(r"[-+]?\d*\.\d+|\d+", t_score) # Extract float values


	return(f_score, t_score)


def df_gen(f_data, t_data, gene_title, f_seqs, t_seqs, df_title):
	""" 
	Dataframe generating.
	Function takes an input 5' and 3' ss scores list (f_data, t_data),
	input GenBank file name (gene_title),
	seq of 5' and 3' ss (f_seqs, t_seqs) and output dataframe name;
	Function returns CSV file

	"""
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

with open('seq_log', 'w') as file:  # Initialiazing log file
		file.write('Used features details and secondary info: \n \n')

with open(df_name, 'w', newline = '') as res_file:  # Initializing total df
	writer = csv.writer(res_file)
	writer.writerow(['score', 'int_num', 'end', 'site_seq', 'gene'])

print('Genes progress:')
for file in glob.glob('*.gb'):
	print(file.split('.')[0])
	test_record = file
	f_ss_sequences, t_ss_sequences = ss_extract(test_record)
	f_ss_scores, t_ss_scores = maxent_calc(f_ss_sequences, t_ss_sequences,
		                                   test_record)
	df_gen(f_ss_scores, t_ss_scores, test_record,
		   f_ss_sequences, t_ss_sequences, df_name)


# That's all!