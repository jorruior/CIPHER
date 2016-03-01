# CIPHER
Online tool: http://evolutionarygenomics.imim.es/josepl/pol_web/index.php/primera/home

CIPHER is a program for the prediction of coding sequences, or open reading frames (ORFs), in transcripts. It is based on the differences in hexamer nucleotide frequencies in sets of known coding and non-coding sequences for different species (see Algorithm). It was used for the first time in the analysis of putative coding transcripts identified by ribosome profiling (Ruiz-Orera et al., 2014).  The user should upload a set of DNA sequences in FASTA format. Possible uses of the program are the prediction of putative proteins in a set of assembled transcripts or the detection of additional peptides produced from a given mRNA. The analysis is restricted to ORFs starting with an ATG codon (Methyonine) and ending with a STOP codon. Empirical p-value cutoffs have been precalculated for ORFs of different minimum length in various model species. 


MAIN EXECUTABLE: cipher.py (-h for options)
Required parameters:
-i/--input: Transcript file in FASTA format.
-o/--outfile: Name of two output files (orfs and table).
-s/--specie: Species: human,primate,mouse,drosophila,arabidopsis,yeast.
-x/--table: Python object with precomputed discores. Use dicodon_matrices.py to build novel models.
Other parameters:
-p/--pvalue: P-value based on a negative model in random introns. 0.05 (default and recommended), 0.01, 0.005
-t/--threshold: Minimum predicted ORF length (amino acids). 24 (default) or 60.
-n/--orf_number: Select in 'longest' (default) or 'all' ORFs per transcript are analysed.

Example to run the test on cipher:
python cipher.py -i tutorial/sequences_human.fa -o test/test -s human -x tables/hsa_coding_to_intron_dicodon_usage.obj         

FILE TO CREATE DICODONS TABLES: dicodon_matrices.py  (-h for options)

scores.py contains a function with computed thresholds for different ORFs and p-values.


References:
Long non-coding RNAs as a source of new peptides Ruiz-Orera, J., Messeguer, X., Subirana, J. A., & Alba, M. M. (2014). Long non-coding RNAs as a source of new peptides. eLife, 3, 1 ‒ 24.



