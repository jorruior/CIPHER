#!/usr/bin/env python
'''---CIPHER---'''

import pickle
import sys
from optparse import OptionParser
import os
##Modules of the script
import scores

__author__ = "Jorge Ruiz-Orera"
__contributor__="Jorge Ruiz-Orera, Pol Verdaguer-Grau, Jose Luis Villanueva-Canyas, Xavier Messeguer, M.Mar Alba"
__copyright__ = ""
__credits__ = []
__license__ = ""
__version__="1.0.0"
__maintainer__ = "Jorge Ruiz-Orera"
__email__ = "jruiz@imim.es; malba@imim.es"


class orf_object:
	def __init__(self, sequence, start, end):
		self.sequence = sequence
		self.start = start
		self.end = end

def find_all(sequence, subsequence):
	''' Returns a list of indexes within sequence that are the start of subsequence'''
	start = 0
	idxs = []
	next_idx = sequence.find(subsequence, start)
 
	while next_idx != -1:
		idxs.append(next_idx)
		start = next_idx + 1     # Move past this on the next time around
		next_idx = sequence.find(subsequence, start)
 
	return idxs
		
def find_orfs(sequence, threshold):
	''' Finds all valid open reading frames and returns them as a list with start and stop relative coordinates'''
 
	starts = find_all(sequence, 'ATG')
	stop_amber = find_all(sequence, 'TAG')
	stop_ochre = find_all(sequence, 'TAA')
	stop_umber = find_all(sequence, 'TGA')
	stops = stop_amber + stop_ochre + stop_umber
	stops.sort()
 
	orfs = []
 
	for start in starts:
		for stop in stops:
			if start < stop \
				and (start - stop) % 3 == 0:
					if len(sequence[start:stop+3]) >= int(threshold)*3: #Only ORFS > threshold
						orf_obj = orf_object(sequence[start:stop+3], start, stop+3)
						orfs.append(orf_obj)
					break

	orfs.sort(key=lambda x: len(x.sequence), reverse=True)
	return orfs

def find_cod_orfs(sequence,name,threshold,ORF_s,dicodons,output,output2,coding_threshold):
	''' Compute scores and write output in predicted ORFs ''' 
	orfs = find_orfs(sequence, threshold)
	ends = []
	counter = 0
	orfs_n = 1
	for orf in orfs:
		cds_seq = orf.sequence
		if (counter > 0) and (ORF_s == "longest"):
			break
		if orf.end in ends:
			continue
		annot_discore = 0
		for (n,i) in enumerate(xrange(3, len(cds_seq)-3, 3)): #Skipping first ATG and last stop
			dicodon = cds_seq[i:i+6]
			try:
				annot_discore = annot_discore + float(dicodons[dicodon])
			except:
				pass

		annot_discore = float(annot_discore)/float(n)
		if ((len(str(cds_seq)) >= 300) & (annot_discore >= coding_threshold[0])) | ((len(str(cds_seq)) < 300) & (len(str(cds_seq)) >= 180) & (annot_discore >= coding_threshold[1])) | (len(str(cds_seq)) < 180) & (annot_discore >= coding_threshold[2])): 
			ends.append(orf.end)
			output.write(name + "\torf_" + str(orfs_n) + "\t" + str(orf.start) + "-" + str(orf.end) + "\t" + str(len(str(cds_seq))) + "\t" + str(len(sequence)) + "\t" + str(format(annot_discore, '.5f')) + "\n")
			output2.write("> " + name + "_" + str(orfs_n) + "\n" + str(cds_seq) + "\n")
			counter = 1
			orfs_n += 1

	return counter

def check_arg (arg_str,s):
	'''Check if arg was written'''
	if not arg_str:   # if filename is not given
		print "Error: " + str(s) + " argument not given\n"
		exit()

def check_file (file_str):
	'''Check if input really exists'''
	try:
		open("%s" %file_str)
	except:
		print "Error: " + file_str + " input not found\n"
		exit()

def main():
	usage = "\n%prog  [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input",action="store",dest="fasta",help="Transcript file in FASTA format (required).")
	parser.add_option("-o","--outfile",action="store",dest="out_file",help="Output file (required). Tab separated text file: name <tab> ORF number <tab> ORF position <tab ORF length <tab> Transcript length <tab> Coding Score.")
	parser.add_option("-s","--specie",action="store",dest="sp",help="Specie to be analyzed (required). Prebuilt thresholds and codon matrices for: human, primate, mouse, zebrafish, arabidopsis, yeast")
	parser.add_option("-x","--table",action="store",dest="dicodons",help="Python object with dicodon scores (required)")
	parser.add_option("-p","--pvalue",action="store",dest="p_value",default="0.05",help="P-value: 0.05, 0.01, 0.005, none. default=0.05")
	parser.add_option("-t","--threshold",action="store",type=int,dest="threshold",default="24",help="Minimum ORF length threshold in amino acids (>4). default=24")
	parser.add_option("-n","--orf_number",action="store",dest="ORF_s",default='longest',help="Predict longest or all ORFs. default=longest")
	
	(opt,args)=parser.parse_args()

	check_arg(opt.fasta,"--input")
	check_arg(opt.out_file,"--outfile")
	check_arg(opt.sp,"--specie")
	check_arg(opt.dicodons,"--table")
	check_file(opt.fasta)
	check_file(opt.dicodons)

	if int(opt.threshold) < 4:
		print("Minimum ORF threshold should be an integer larger than 4")
		exit(0)
		
	(specie,coding_threshold) = scores.find_scores(opt.sp,opt.p_value)

	dicodons = pickle.load(open(opt.dicodons,'r'))

	output = open(opt.out_file + "_scores.txt","w+")
	output.write("sequence\tORF_number\tORF_pos\tORF_len\ttranscript_len\tcoding_score\n")
	output2 = open(opt.out_file + "_orfs.fa","w+")

	c = 0
	t = 0
	with open(opt.fasta) as fp:
		for name, seq in scores.read_fasta(fp):
			c = c + find_cod_orfs(seq.upper(),name.replace(">","").replace(" ",""),opt.threshold,opt.ORF_s,dicodons,output,output2,coding_threshold)
			t += 1

	print >>sys.stderr, "Sequences with coding potential: " + str(c) + " out of " + str(t) + " (" + str(round(float(c)/float(t)*100,2)) + "%)"
	output.close()
	output2.close()


if __name__ == '__main__':
	main()

exit(0)
