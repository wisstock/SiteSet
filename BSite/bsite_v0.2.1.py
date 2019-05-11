#!/usr/bin/env python3

import sys
import csv
import glob
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import Align


"""
Localizing of semi-specific binding sites (bs) of RNA-binding proteins (RBP)
using local alignment algorithm

Copyright © 2018 Borys Olifirov
Institure of Molecular Biology and Genetics of NASU, 
Systems Biology Research Group

Script return CSV dataframe bsites_df.csv
('gene', 'factor', 'site', 'score', 'location')
and log file bsite.log

The script require Python 3.5 or higer and Biopyton 1.68 or higer


v0.1 - search semi-seq specific binding sites and generating data frame
v0.2 - recognize target intron and localize alignment only in their seq
v0.2.1 - add int_start and int_end columns to output data frame
"""


def amb_nuc(seq):
	"""
	Interprets ambiguous nucleotide symbols 
	and generating all possible sequences 
	for all nucleotide positional variations

	""" 
	combinations = []
	ambiguity_dict = {
	    'A':['A'],
	    'C':['C'],
	    'G':['G'],
	    'T':['T'],
	    'W':['A','T'],
	    'S':['C','G'],
	    'M':['A','C'],
	    'K':['G','T'],
	    'R':['A','G'], 
	    'Y':['C','T'],
	    'B':['C','G','T'],
	    'D':['A','G','T'],
	    'H':['A','C','T'],
	    'V':['A','C','G'],
	    'N':['A','C','G','T'],
	    'Z':[]
	}

	# W     Weak
	# S     Strong
	# M     aMino
	# K     Keto
	# R     puRine
	# Y     pYrimidine
	# B     not A (B comes after A)
	# D     not C (D comes after C)
	# H     not G (H comes after G)
	# V     not T (V comes after T and U)
	# N     any Nucleotide (not a gap)
	# Z     Zero

	for symbol in seq:
	    if symbol == 'Z':
	        print('Zero symbol encountered. Skipping.    ')
	    new_combinations = []
	    if combinations:    
	        for combination in combinations:
	            for nucleotide in ambiguity_dict[symbol]:
	                new_combinations.append(combination + nucleotide)
	    else:
	        for nucleotide in ambiguity_dict[symbol]:
	            new_combinations.append(nucleotide)
	    if new_combinations:
	        combinations = new_combinations


	return combinations


def intron_ext(gene_file, apa_location):
	"""
	Selective exrtraction of introns with experemental proved polyadenelation

	Function takes an input GB-file name (gene_file) and number
	of target intron as integer (apa_location);
	
	Function returns target intron seq as string

	"""
	i = 1

	seq_record = SeqIO.read(gene_file, 'genbank')
	gene_seq = str(seq_record.seq)

	for feature in seq_record.features:  # Loop over the gen file, 
		if feature.type == 'mRNA':       # get the intron with APA site
			for exon_location in feature.location.parts:
				if i <= apa_location: 
					ex_end = int(exon_location.end)
					i += 1
				else: 
					break

	i = 1

	for feature in seq_record.features:  # Loop over the gen file, 
		if feature.type == 'mRNA':       # get the introne with APA site
			for exon_location in feature.location.parts:
				if i-1 <= apa_location: 
					ex_start = int(exon_location.start)
					i += 1
				else: 
					break

	with open('bsite.log', 'a') as log_file:        # adding data to log file
			log_file.write('Selected feature ID:' + \
				           str(feature.qualifiers.get("db_xref")) + '\n')
			log_file.write('Target intron lenght:' + \
				           str(len(seq_record.seq[ex_end:ex_start])) + '\n')          


	return(seq_record.seq[ex_end:ex_start], ex_end, ex_start)


def algn_fuck(gene_n, ref_seq, start, ends, bss_list, factor, 
	          table_name, aligner_parameters = [2, -5, -0.1, -0.5]):
	"""
	Search and localizing individual target sites using PairwiseAligner;

	Function takes an input reference seq in string format (ref_seq),
	list of binding site variation from amb_nuc function (bss_list) 
	and aligner parameters (aligner_prmtrs) as a vector of score values 
	[match score, mismatch score, open gap score and extend gap score];

	Function return list of aligner parameters scores (scores), 
	list of target's sites first nucleotide locations (starts),
	alignment iterations (algn_it) and aligner parameters (aligner)

	"""
	scores = []
	starts = []
	ind_score = []
	ind_start = []
	algn_it = 0
	start_num = 0

	
	aligner = Align.PairwiseAligner()
	aligner.mode = 'local'
	aligner.target_internal_open_gap_score = -1000
	aligner.target_internal_extend_gap_score = -10000
	aligner.target_left_gap_score = -1000
	aligner.targen_right_gap_score = -1000
	aligner.match, aligner.mismatch, aligner.query_internal_open_gap_score, \
	    aligner.query_internal_extend_gap_score = aligner_parameters


	for bs_seq in bss_list:
		local_key = bs_seq
		alignments = aligner.align(ref_seq, bs_seq)
		alg_list = list(alignments)
		
		for alg in alignments:  # present individual aligment & compute score
			score_num = alg.score
			ind_score.append(score_num)
			str_alg = str(alg)
			str_len = int(len(str_alg)/3)
			str_str = str_alg[str_len*2:str_len*3]
			print('=', sep = '',           # dynamic print progress bar, 
			      end = '', flush = True)  # one '=' per alignment result

			while str_str[start_num] != bs_seq[0]:  # compute site position
				start_num += 1

			one_row = [gene_n, factor, bs_seq, score_num, start, ends, start_num+1]
			with open(table_name, 'a', newline = '') as res_file:
				writer = csv.writer(res_file)
				writer.writerow(one_row)

			ind_start.append(start_num+1)
			start_num = 0

			algn_it += 1


		scores.append(ind_score)
		starts.append(ind_start)

		ind_score = []
		ind_start = []

	with open('bsite.log', 'a') as log_file:        # adding data to log file
		log_file.write('Alignment count:' + str(algn_it) + '\n')
		log_file.write('========================================\n')  # sepp line


	return(aligner)
	


df_name = 'bsite_df.csv'
algnr_prmtrs = [2, -0.1, -1.9, -100]
apa_int = {'GART': 11, 'CSTF3': 3, 'NAP1L': 13, 'PCIF1': 2, 'ZMYM3': 7}


with open('bsite.log', 'w') as log_file:  # initialiazing log file
		log_file.write('========== BSite v0.2 | RESULTS REVIEW ==========\n')
		log_file.write('=================================================\n')

with open(df_name, 'w', newline = '') as res_file:  # initializing total dataframe
	writer = csv.writer(res_file)
	writer.writerow(['gene', 'b_factor', 'site', 'score', 'int_start', 'int_end', 'location'])


print('\n\n========== BSite v0.2 | RESULTS REVIEW ==========')
print('=================================================\n\n')
with open('factors.config', 'r', newline = '') as factors_list:  # read line by line list of factors
	reader = csv.reader(factors_list)
	for row in reader:
		ever_seq = amb_nuc(row[1])
		factor_name = row[0]
		print('===', factor_name, 'BS searching')

		with open('bsite.log', 'a') as log_file:                      # adding data to log file
				log_file.write('\n\n======== Factor name:' + factor_name + '\n\n')           # factor name
				log_file.write('Site`s variants:' + str(len(ever_seq)) + '\n')     # total count of binding site variations
				log_file.write('Seq:' + str(ever_seq) + '\n')                      # sequence for each site's variant

		for gb_file in glob.glob('*.gb'):  # loop over all gb-file, align actual binding site to each target intron
			gene_name = gb_file.split('.')[0]
			int_int = apa_int[gene_name]  # get target introne number by gb-file name
			print('\n', gene_name, 'in progress')

			with open('bsite.log', 'a') as log_file:                  # adding data to log file
				log_file.write('\nGene name:' + gene_name + '\n')
				log_file.write('Target intron num:' + str(int_int) + '\n')

			int_seq, int_start, int_end  = intron_ext(gb_file, int_int)
			aligner_scores = algn_fuck(gene_name, int_seq, int_start, int_end, ever_seq, factor_name, df_name)
			print('\n', gb_file.split('.')[0], 'done\n')

# That's all!