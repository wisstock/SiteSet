#!/usr/bin/env python3

import csv
import numpy as nm
import pandas as pd
import glob
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import Align

"""
Parsing localization of elected binding sites;

Copyright © Borys Olifirov 2018
Institure of Molecular Biology and Genetics of NASU,
Systems Biology Research Group

The script require Python 3.5 or higer and Biopyton 1.68 or higer

"""


class SiteSet(): 
	"""
	Binding site associates with gene sequence and features location
	Method's reports are beginning by '|' character
	('| Factors list success uploated')
 
	"""

	def __init__(self, indep_input, dep_input, factors_file=None, dep_file=None, apa_loc=None, gb_file=None):
		"""
		Upload gene factors lists
		and optional factors list in CSV format, results of semisite_v1.0.py,
		columns name - gene,b_factor,site,score,location

		"""
		self.RelLoc = pd.DataFrame(columns=['gene', 'indep_factor', 'indep_site_position', 'dep_factor', 'Lr'])
		self.Z_sore = pd.DataFrame(columns=['gene', 'indep_factor', 'indep_site_position', 'dep_factor', 'dep_factor_position', 'z_score'])

		self.ApaLoc = apa_loc

		try:
			self.DepList = dep_input
			self.IndepList = indep_input
			self.GeneData = {}
			print('|SiteSet| Facors list was succesfull uploated\n')
		except:
			print('|SiteSet|!| Factors list dosen\'t define! SiteSet instancw wasn\'t created\n')


		try:
			self.FactorsTable = pd.read_csv(factors_file)  # upload total dataframe (SemiSite results)
			print('|SiteSet| SemiSite data frame was successful uploated\n')
		except:
			print('|SiteSet|!| SemiSite data frame donesn\'t define, no input file! Please, upload data later\n')

			
		try:
		    first_gene = SeqIO.read(gb_file, 'genbank')
		    gene_name = gb_file.split('.')[0]
		    self.GeneData = dict([(gene_name, first_gene)])
		    print('|SiteSet| Gene list was was successful initiated,', gene_name, 'GB-file was uploated\n')
		except:
			print('|SiteSet|!| No one GB-files wasn\'t found, you can upload GB-files later\n')


		if len([col for col in list(self.FactorsTable.columns) if col == 'dep']) == 0:  # the presence of dependences (dep) column checking
			indep_list = []
			self.FactorsTable['dep'] = 'dep'  # creating blank "dep" column

			with open(dep_file, 'r', newline = '') as factors_dep:  # reading factors_dep file
				reader = csv.reader(factors_dep)
				for row in reader:
					if row[1] == 'False':
						indep_list.append(row[0])
					else:
						continue
			self.FactorsTable.loc[self.FactorsTable['factor'].isin(self.IndepList), 'dep'] = 'indep'  # replace 'dep' values for independences factors
			print('|SiteSet| Dependencies column was successful generated and uploaded\n')
		else:
			print('|SiteSet| Dependencies values is already avaliable\n')  # DONE
			

	def GetTable(self, write_option=None):
		"""
		Write FactorsTable with dep column as CSV file

		"""
		
		if write_option == 'factors':
			pass
		elif write_option == 'z-score':
			pass
		elif write_option == 'relloc':
			pass
		elif write_option == 'all':
			pass
		else:
			pass

		print(self.FactorsTable)
		print(self.FactorsTable.loc[self.FactorsTable['factor'].isin(['PAS', 'SRSF3']), ['factor', 'dep']])
		# print('Rows count ' + str(len(self.FactorsTable.index)))  # WHY?


	def AddGene(self, new_file):
		"""
		Appending GB-file to gene dictionary, 
		GB-file name will use as a key and SeqIO result as a value

		"""
		if len([a for a in list(self.GeneData.keys()) if a == new_file.split('.')[0]]) == 0:  # gene duplication checking
		    self.GeneData.update([(new_file.split('.')[0], SeqIO.read(new_file, 'genbank'))])
		    print('|AddGene| Gene list was was successful initiated,', new_file.split('.')[0], 'GB-file was uploated\n')
		else:
			print('|AddGene|!|', new_file.split('.')[0], 'is already available in GeneData')  # DONE


	def GetExon(self, gene=None, option='count'):
		"""
		Depending on the option parameter method returns exons count 
		(default value, 'count'),
		first exon seq ('first'), last exon seq ('last') or exon seq by number

		"""

		ex_f_loc = []
		ex_t_loc = []
		
		# Loop over the gen file, get the mRNA starts and ends position
		# for '+' strand, detecting first mRNA-type feature
		# and extract 5'-ends (f_loc) and 3'-ends locations of annotated exon
		for feature in self.__row_gen_record.features:
			if feature.type == 'mRNA':
				for exon_location in feature.location.parts:
					ex_f_loc.append(int(exon_location.start))
					ex_t_loc.append(int(exon_location.end))
				break

		if option == 'count':
			print('| ' + self.GeneName + ' counts of ' + str(len(f_loc)) + \
				  ' exons\n| Enter exon number and try again\n')
		elif option == 'first':
			return self.GeneSeq[ex_f_loc[0] : ex_t_loc[0]]
		elif option == 'last':
			return self.GeneSeq[ex_f_loc[len(ex_f_loc) - 1] : ex_t_loc[len(ex_f_loc) - 1]]
		elif type(option) == int:
			return self.GeneSeq[ex_f_loc[option - 1] : ex_t_loc[option - 1]]


	def GetIntron(self, gene=None, option = 'count'):
		"""
		Depending on the option parameter method returns introns count 
		(default value, 'count'),
		first intron seq ('first'), last intron seq ('last') 
		or intron seq by number

		"""
		
		if len([a for a in list(self.GeneData.keys()) if a == gene]) == 0:  # gene checking
			print('|GetIntron|!| No match found, try another gene')
		else:
			__local_gene = self.GeneData[gene]

			in_f_loc = []
			in_t_loc = []
			
			# Loop over the gen file, get the mRNA starts and ends position
			# for '+' strand, detecting first mRNA-type feature
			# and extract 5'-ends (f_loc) and 3'-ends locations of annotated exon
			for feature in __local_gene.features:
				if feature.type == 'mRNA':
					for exon_location in feature.location.parts:
						in_t_loc.append(int(exon_location.start))
						in_f_loc.append(int(exon_location.end))
					break

			if option == 'count':
				print('|GenIntron|', gene, 'counts of ', str(len(in_f_loc)),'introns\n')

			elif option == 'first':
				return self.GeneSeq[in_f_loc[0] : in_t_loc[1]]

			elif option == 'last':
				return self.GeneSeq[in_f_loc[len(in_f_loc) - 2] : in_t_loc[len(in_f_loc) - 1]] 

			elif type(option) == int:
				self.TargetIntStart = in_f_loc[option - 1]
				self.TargetIntEnd = in_t_loc[option]
				return __local_gene.seq[in_f_loc[option - 1] : in_t_loc[option]]  # DONE


	def RelativeLoc(self, feature_type=None, write_output=False):
		"""
		Calculate relative locations of depedent factor sites
		about independent factor
		ans Z-score for each dep sites;
		the search zone may be limited
		to a certain feature (see GetIntron/Exon methods);
		feature_type get an inpup exon/intron value
		and feature_number variable stores the number of the feature.

		"""

		# def erf(x):
		# 	"""
		# 	Error function

		# 	"""

		#     # save the sign of x
		#     sign = 1 if x >= 0 else -1
		#     x = abs(x)

		#     # constants
		#     a1 =  0.254829592
		#     a2 = -0.284496736
		#     a3 =  1.421413741
		#     a4 = -1.453152027
		#     a5 =  1.061405429
		#     p  =  0.3275911

		#     # A&S formula 7.1.26
		#     t = 1.0/(1.0 + p*x)
		#     y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x)
		#     return sign*y # erf(-x) = -erf(x)



		def __relloc_calc(input_indep, input_dep, indep_name, dep_name):
			"""
			Relative location calculating; 
			Fun get as input two lists, locations 
			of dependencies and dependencies binding sites
             
              2  sum(d-i)^2
			Lr = -----------
                    n

                 d-i
            z = -----
                 Lr

			"""

			__lr_res_list, __z_res_list = [], [] 
			__lr_columns_name = ['gene', 'indep_factor', 'indep_site_position', 'dep_factor', 'Lr']
			__z_columns_name = ['gene', 'indep_factor', 'indep_site_position', 'dep_factor', 'dep_factor_position', 'z_score']

			__lr_calc = lambda x: pow(sum([a ** 2 for a in x])/len(x), .5)  # Lr calc function

			for indep_site in input_indep:

				if indep_name == 'PAS' and dep_name == 'up_PAS':  # one-sided calculation for upstream motif sites
					up_dep_input = [abs(x - indep_site) for x in input_dep if x > indep_site]
					try:
						lr = __lr_calc(up_dep_input)
					except ZeroDivisionError:
						print('| No upstream PAS site found!')
						__lr_res_list.append([self.current_gene, indep_name, indep_site, dep_name, None])
						continue
					__lr_res_list.append([self.current_gene, indep_name, indep_site, dep_name, lr])  # add current row to 2D list
					for dep_site in input_dep:  # loop over all dep sites for z-score calculation
						zz = (dep_site - indep_site)/lr
						__z_res_list.append([self.current_gene, indep_name, indep_site, dep_name, dep_site, zz])
				elif indep_name == 'PAS' and dep_name == 'down_PAS':  # one-sided calculation for downstream motif sites
					down_dep_input = [abs(x - indep_site) for x in input_dep if x < indep_site]
					try:
						lr = __lr_calc(down_dep_input)
					except ZeroDivisionError:
						print('| No upstream PAS site found!')
						__lr_res_list.append([self.current_gene, indep_name, indep_site, dep_name, None])
						continue
					__lr_res_list.append([self.current_gene, indep_name, indep_site, dep_name, lr])  # add current row to 2D list
					for dep_site in input_dep:  # loop over all dep sites for z-score calculation
						zz = (dep_site - indep_site)/lr
						__z_res_list.append([self.current_gene, indep_name, indep_site, dep_name, dep_site, zz])
				else:
					normal_dep_input = [abs(x - indep_site) for x in input_dep]  # loc normalization by indep site position
					try:
						lr = __lr_calc(normal_dep_input)
					except ZeroDivisionError:
						print('| No dep site found!')
						__lr_res_list.append([self.current_gene, indep_name, indep_site, dep_name, None])
						continue
					__lr_res_list.append([self.current_gene, indep_name, indep_site, dep_name, lr])  # add current row to 2D list

					for dep_site in input_dep:  # loop over all dep sites for z-score calculation
						zz = (dep_site - indep_site)/lr
						__z_res_list.append([self.current_gene, indep_name, indep_site, dep_name, dep_site, zz])
					a = __z_res_list
						

			return(pd.DataFrame(__lr_res_list, columns=__lr_columns_name), pd.DataFrame(__z_res_list, columns=__z_columns_name))



		if feature_type == None:  # extract indep factor sites location for selected gene
			for in_fac in indep_factor:
				self.__inter_set = self.FactorsTable.loc[self.FactorsTable['factor'] == in_fac]
				self.__in_fac_locs = self.__inter_set.loc[self.__inter_set['gene'] == self.GeneName, 'location']

				for dep_fac in dep_factor:  # extract dep factor sites location for selected gene
					self.__depter_set = self.FactorsTable.loc[self.FactorsTable['factor'] == dep_fac]
					self.__dep_fac_locs = list(self.__depter_set.loc[self.__depter_set['gene'] == self.GeneName, 'location'])

			if write_dataframe == True:
				__output_lr_name = 'RelLoc_' + self.GeneName + '_full.csv'
				self.RelLoc_instance.to_csv(__output_lr_name)
			else:
				print('\n| RelativeLoc CSV data frame wasn\'t created (write_dataframe = False)')


		elif feature_type == 'intron':
			for self.current_gene in list(self.GeneData.keys()):

				int_num = self.ApaLoc[self.current_gene]
				print(list(self.GeneData.keys()))
				print(self.current_gene)

				self.GetIntron(self.current_gene, int_num)  # using GetIntron method for extraction location of target intron ends
				print('|RelLoc|', int_num, 'th intron\'s length is', str(self.TargetIntEnd - self.TargetIntStart), '\n')
				print('|RelLoc| Intron edges', self.TargetIntStart, self.TargetIntEnd, '\n')

				for in_fac in self.IndepList:  # extract indep factor sites location for selected gene forn FactorsTable
					self.__inter_set = self.FactorsTable.loc[self.FactorsTable['factor'] == in_fac]                         #
					self.__in_fac_locs = list(self.__inter_set.loc[self.__inter_set['gene'] == self.current_gene, 'location'])  # too long
					target_indep_sites = [a for a in self.__in_fac_locs if self.TargetIntStart <= a <= self.TargetIntEnd]   # construction
					print('\n|||', in_fac, 'indep sites positions', target_indep_sites)
					

					for dep_fac in self.DepList:  # extract dep factor sites location for selected gene forn FactorsTable
						if dep_fac == in_fac:  # exclude duplicates if same factors counts in dep and indep lists simultaneously
							continue
						else:
							self.__depter_set = self.FactorsTable.loc[self.FactorsTable['factor'] == dep_fac]
							self.__dep_fac_locs = self.__depter_set.loc[self.__depter_set['gene'] == self.current_gene, 'location']
							target_dep_sites = [a for a in self.__dep_fac_locs if self.TargetIntStart <= a <= self.TargetIntEnd]
							print('\n|', dep_fac, 'dep sites positions', target_dep_sites)

							lr_instance, z_instance = __relloc_calc(target_indep_sites, target_dep_sites, in_fac, dep_fac)
							self.RelLoc = self.RelLoc.append(lr_instance, ignore_index=True)  # add results of each dep fac to total DF
							self.Z_sore = self.Z_sore.append(z_instance, ignore_index=True)


			if write_output == True:
				self.RelLoc.to_csv('RelLoc.csv')

				self.Z_sore.to_csv('ZLoc.csv')
			else:
				print('| RelativeLoc CSV data frame wasn\'t created (write_output = False)\n')


		elif feature_type == 'exon':
			pass


	def MinimalLr(self):
		"""
		searching dep site with minimal Lr sum 

		"""

		min_sum, res_sum = [], [] 
		res_pos = []
		min_indep_pos = {}

		for self.current_gene in list(self.GeneData.keys()):
				int_num = self.ApaLoc[self.current_gene]
				self.GetIntron(self.current_gene, int_num)

				for indep_fac in self.IndepList:
					instance = self.FactorsTable.loc[(self.FactorsTable['factor'] == indep_fac) & \
					    (self.FactorsTable['gene'] == self.current_gene), \
					    'location']

					for i in [a for a in instance if self.TargetIntStart <= a <= self.TargetIntEnd]:

						res_sum = self.RelLoc.loc[(self.RelLoc['indep_site_position'] == i) & (self.RelLoc['gene'] == self.current_gene),'Lr']
						res_sum = sum([a for a in res_sum if a !=None])
						min_sum.append(res_sum)
						res_pos.append(i)
					print(self.current_gene, res_pos[min_sum.index(min(min_sum))])
				min_sum, res_pos = [], []


	def IndepCompare(self, compare_name, dep_factor, indep_factor):
		"""
		Comparing relative locations of dependent factor between two genes

		"""
		pass


	def ApaDetector(self, canonic = True):
		"""
		Detecting alternative polyadenelation-related motifs
		in correct mutual arrangement
		== 5'ss ==> up_PAS ==> PAS ==> down_PAS ==> 3'ss ==

		If 'canonic' True search will limit only introns seq 
		(False - full gen searching)

		"""
		pass
		


dep_list = ['SRSF3', 'SRSF10']
indep_list = ['YTHDC1']

apa = apa_int = {'GART': 11, 'CSTF3': 3, 'NAP1L': 13, 'PCIF1': 2, 'ZMYM3': 7}

m6a = SiteSet(indep_list, dep_list, 'semisite_df.csv', 'factor_dep.config', apa)

for file in glob.glob('*.gb'):
	m6a.AddGene(file)

# print(m6a.FactorsTable)


m6a.RelativeLoc('intron', True)
m6a.NextDoor()
			
	

# That's all!