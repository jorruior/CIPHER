#!/usr/bin/env python
'''---Build dicodon matrices for CIPHER---'''

import sys
import os
import math
import pickle
from optparse import OptionParser
##Modules of the script
import scores

__author__ = "Jorge Ruiz-Orera"
__contributor__="Jorge Ruiz-Orera, Pol Verdaguer-Grau, Jose Luis Villanueva-Canyas, Xavier Messeguer, M.Mar Alba"

def do_dict():
	d = {}
	f = {}
	nucleotides = ("A","C","G","T")
	for nt in nucleotides:
		for nt2 in nucleotides:
			for nt3 in nucleotides:
				for nt4 in nucleotides:
					for nt5 in nucleotides:
						for nt6 in nucleotides:
							d[(nt+nt2+nt3+nt4+nt5+nt6,"c")] = 0
							d[(nt+nt2+nt3+nt4+nt5+nt6,"nc")] = 0
							f[nt+nt2+nt3+nt4+nt5+nt6] = 0							
	return d,f

def score_orfs(cds_seq,name,t):
	total_codons = 0
	for (n,i) in enumerate(xrange(3, len(cds_seq)-3, 3)): #Skipping first ATG and last stop
		dicodon = cds_seq[i:i+6]
		try:
			codons[(dicodon,t)] += 1
			total_codons += 1
		except:
			pass
	return total_codons

def main():
	usage = "\n%prog  [options]"
	parser = OptionParser(usage)
	parser.add_option("-c","--coding",action="store",dest="c",help="Coding file: ORFs in nucleotide FASTA format (required).")
	parser.add_option("-n","--noncoding",action="store",dest="nc",help="Non-coding file: ORFs in nucleotide FASTA format (required).")
	parser.add_option("-o","--outfile",action="store",dest="out_file",help="Output file, written in 'tables' directory (required).")

	(opt,args)=parser.parse_args()

	global codons
	(codons,final_codons) = do_dict()

	c = 0
	with open(opt.c) as fp:
		for name, seq in scores.read_fasta(fp):
			c = c + score_orfs(seq,name.replace(">","").replace(" ",""),"c")

	nc = 0
	with open(opt.nc) as fp:
		for name, seq in scores.read_fasta(fp):
			nc = nc + score_orfs(seq,name.replace(">","").replace(" ",""),"nc")

	for key in final_codons.keys():
		if float(codons[(key,"c")] > 0) and float(codons[(key,"nc")] > 0):
			final_codons[key] = math.log10((float(codons[(key,"c")]) / float(c))/(float(codons[(key,"nc")]) / float(nc)))
		elif float(codons[key] > 0): #pseudocount
			final_codons[key] = math.log10((float(codons[(key,"c")]) / float(c))/(1 / float(nc)))
		elif float(codons_i[key] > 0): #pseudocount
			final_codons[key] = math.log10((1 / float(c))/(float(codons[(key,"nc")]) / float(nc)))
		else:
			final_codons[key] = 0

	file_obj = open("./tables/" + opt.out_file, 'w')
	pickle.dump(final_codons, file_obj)


if __name__ == '__main__':
	main()

exit(0)
