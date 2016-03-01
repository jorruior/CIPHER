#!/usr/bin/env python

__author__ = "Jorge Ruiz-Orera"
__contributor__="Jorge Ruiz-Orera, Pol Verdaguer-Grau, Jose Luis Villanueva-Canyas, Xavier Messeguer, M.Mar Alba"

def find_scores(sp,p_value):
	''' Finds thresholds based on a negative model in intron random ORFs '''
	if (sp == "human") or (sp == "primate") or (sp == "hsa"):
		specie = "hsa"
		if p_value == "0.05":
			coding_threshold = (0.0448012743458, 0.0313544171427, 0.022090633722)       
		elif p_value == "0.01":
			coding_threshold = (0.103283772197, 0.0828501857409, 0.0788118157357)
		elif p_value == "0.005":
			coding_threshold = (0.128866695709, 0.105405369318, 0.107772536839)
		  
	elif (sp == "mouse") or (sp == "mmu"):
		specie = "mmu"
		if p_value == "0.05":
			coding_threshold = (0.0700036382617, 0.0511031277081, 0.0346298615614)
		elif p_value == "0.01":
			coding_threshold = (0.138588337003, 0.109739472571, 0.101682499377)           
		elif p_value == "0.005":
			coding_threshold = (0.169838369021, 0.139522454066, 0.151992160178)

	elif (sp == "zebrafish") or (sp == "dan"):
		specie = "dan"  
		if p_value == "0.05":
			coding_threshold = (0.04765878026, 0.0337564755358, 0.030315191054)
		elif p_value == "0.01":
			coding_threshold = (0.138144489602, 0.104594395726, 0.11826576561)			
		elif p_value == "0.005":
			coding_threshold = (0.17579879022, 0.137402931564, 0.151342956871)	
		  
	elif (sp == "arabidopsis") or (sp == "ath"):
		specie = "ath"
		if p_value == "0.05":
			coding_threshold = (0.0197578288187, -0.0169329278932, -0.02526262537)
		elif p_value == "0.01":
			coding_threshold = (0.0956191042192, 0.0653912502108, 0.0998476799308)
		elif p_value == "0.005":
			coding_threshold = (0.139637050944, 0.120397034526, 0.184112277183)
		  
	elif (sp == "yeast") or (sp == "sac"):
		specie = "sac"
		if p_value == "0.05":
			coding_threshold = (0.1013131163, 0.0562096181754, -0.00291471505818)
		elif p_value == "0.01":
			coding_threshold = (0.171938901579, 0.114368235407, 0.0427984903964)
		elif p_value == "0.005":
			coding_threshold = (0.205835816424, 0.151408870679, 0.0658027186213)

	return specie,coding_threshold

def read_fasta(fp):
	''' Parse a Fasta file, based on Biopython parser'''
	name, seq = None, []
	for line in fp:
		line = line.rstrip()
		if line.startswith(">"):
			if name: yield (name, ''.join(seq))
			name, seq = line, []
		else:
			seq.append(line)
	if name: yield (name, ''.join(seq))
