#!/usr/bin/env python3

"""
Copyright © 2018-2019 Institure of Molecular Biology and Genetics of NASU, Systems Biology Research Group
Copyright © 2018-2019 Borys Olifirov

Localizing of semi-specific binding sites (bs) of RNA-binding proteins (RBP)
using local alignment algorithm.

Script take an input GenBank files with genes 
and CSV file "factors.config" without title row;
with list of factors and their binding sites.

Script return CSV dataframe bsites_df.csv 
('gene', 'factor', 'site', 'score', 'location').

The script require Python 3.5 or higer and Biopyton 1.68 or higer

Consensus binding sites for selected RBPs:
YTHDC1		GGACH
SRSF3		CHWCHMC
SRSF10		TVAAGAHY

v1.1.0.1 - added loging
v1.1.0.2 - intron seq extraction (modified BSite function)

"""


import sys
import csv
import logging
import glob

import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import Align



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


def get_intron(gene=None, option='count'):
	"""
	Depending on the option parameter method returns introns count 
	(default value, 'count'),
	first intron seq ('first'), last intron seq ('last') 
	or intron seq by number

	"""
	
	try:
		seq_record = SeqIO.read(gene, 'genbank')
		logging.debug('seq for intron extraction uploated')
	except:
		logging.error('gene not define!')

	try:
		type(option) == int
	except:
		logging.error('intron number not define!')

	in_f_loc = []
	in_t_loc = []
	i = 1
		
	# Loop over the gen file, get the mRNA starts and ends position
	# for '+' strand, detecting first mRNA-type feature
	# and extract 5'-ends (f_loc) and 3'-ends locations of annotated exon
	for feature in seq_record.features:
		if feature.type == 'mRNA':
			logging.debug('feature select '+str(feature.type))
			logging.debug('feature select '+str(feature.qualifiers.get("transcript_id")))
			for exon_location in feature.location.parts:
				if i <= option + 1:
					in_t_loc.append(int(exon_location.start))
					in_f_loc.append(int(exon_location.end))
					i += 1
					logging.debug('feature in work '+str(feature.type)+' '+str(feature.qualifiers.get("transcript_id")))
				else:
					break

	target_int_seq = seq_record.seq[in_f_loc[option-1] : in_t_loc[option]]

	logging.info('intron '+str(i-2)+' seq extracted succesfull')
	logging.info('selected feature ID:'+str(feature.qualifiers.get("transcript_id")))
	logging.info('intron lenght: '+str(len(target_int_seq)))

	return target_int_seq



def algn_fuck(gene_n=None,
		      factor=None,
		      seq_mode=None,
		      int_seq=None,
	          ref_seq=None,
	          bss_list=None,
	          table_name=None,
	          aligner_parameters=[2, -5, -0.1, -0.5],
	          ):
	"""
	Search and localizing individual target sites using PairwiseAligner;
	Function takes an input reference seq in string format (ref_seq), list of
	binding site variation from amb_nuc function (bss_list) and aligner
	parameters (aligner_prmtrs) as a vector of score values 
	[match score, mismatch score, open gap score and extend gap score];
	Function return list of aligner parameters scores (scores),
	list of target's sites first nucleotide locations (starts),
	alignment iterations (algn_it) and aligner parameters (aligner)

	"""


	if type(seq_mode) != int:
		seq_record = SeqIO.read(ref_seq, 'genbank')
		gene_seq = str(seq_record.seq)
		logging.debug(len(gene_seq))
	else:
		gene_seq = int_seq

		


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
		alignments = aligner.align(gene_seq, bs_seq)
		alg_list = list(alignments)
		
		for alg in alignments:  # present individual aligment and compute score
			score_num = alg.score
			ind_score.append(score_num)
			str_alg = str(alg)
			str_len = int(len(str_alg)/3)
			str_str = str_alg[str_len*2:str_len*3]

			while str_str[start_num] != bs_seq[0]:  # compute site position
				start_num += 1

			one_row = [gene_n, factor, bs_seq, score_num, start_num+1]
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

	return(scores, starts, algn_it, aligner)	
	



logging.basicConfig(level=logging.DEBUG, filename = 'semisite.log', filemode = 'w')
	# format = '%(levelname)-8s [%(asctime)s] %(message)s', level = logging.DEBUG, filename = 'semisite_%(asctime)s.log', filemode= 'w'
	# debug
	# info
	# warning
	# error
	# critical


df_name = 'semisite_df.csv'                                              #
algnr_prmtrs = [2, -0.1, -1.9, -100]                                     # input data
apa_int = {'GART': 11, 'CSTF3': 3, 'NAP1L': 13, 'PCIF1': 2, 'ZMYM3': 7}  #


with open(df_name, 'w', newline='') as res_file:  # initializing
	writer = csv.writer(res_file)                 # total dataframe
	writer.writerow(['gene', 'factor', 'site', 'score', 'location'])


print('\n\n========== SemiSite v1.0 progress ==========')
print('============================================\n\n')

with open('factors.config', 'r', newline='') as factors_list:  # cinfig file reading
	reader = csv.reader(factors_list)

	for row in reader:
		ever_seq = amb_nuc(row[1])
		factor_name = row[0]

		print('\n\n=========', factor_name, 'BS searching\n')

		logging.info(factor_name + ' BS searching')                   # factor in progress
		logging.info('binding site variants: '+str(ever_seq)+'\n\n')           # print all site variants


		for gb_file in glob.glob('*.gb'):
			gene_name = gb_file.split('.')[0]
			int_int = apa_int[gene_name]  # select current gene target intron namber

			logging.info(gene_name + ' start\n')                         # gene start

			print('===', gene_name, 'in progress')

			get_intron(gb_file, int_int)
			score_list, start_list, algn_num, aligner_scores = algn_fuck(gene_n = gene_name,
				                                                         factor = factor_name,
				                                                         ref_seq = gb_file,
				                                                         bss_list = ever_seq,
				                                                         table_name = df_name
				                                                         )


			print('===', gb_file.split('.')[0], 'done\n')

			logging.info('alignment count: '+str(algn_num)+'\n')             # alignment count
			logging.info(gene_name + ' done\n')                         # gene done
		logging.info('aligner settings:\n'+str(aligner_scores))     # alignment settings

# That's all!