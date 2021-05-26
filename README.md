# CIPHER

CIPHER is a python executable for the prediction of coding sequences, or open reading frames (ORFs), in transcripts. It is based on the differences in hexamer nucleotide frequencies in sets of known coding and non-coding sequences for different species (see Algorithm in the webpage version). It was used for the first time in the analysis of putative coding transcripts identified by ribosome profiling (Ruiz-Orera et al., 2014).  The user should upload a set of DNA sequences in FASTA format. Possible uses of the program are the prediction of putative proteins in a set of assembled transcripts or the detection of additional peptides produced from a given mRNA. The analysis is restricted to ORFs starting with an ATG codon (Methyonine) and ending with a STOP codon. Empirical p-value cutoffs have been precalculated for ORFs of different minimum length in various model species. 


***MAIN EXECUTABLE:*** cipher.py (-h for help)
Required parameters:
```
-i/--input: Transcript file in FASTA format.

-o/--outfile: Name of two output files: orfs and table.

-s/--specie: Species: human,primate,mouse,drosophila,arabidopsis,yeast.

-x/--table: Python object with precomputed discores. Already prebuilt tables available in 'tables' directory for human,mouse,drosophila,arabidopsis and yeast. Use dicodon_matrices.py (-h for help) to build novel models.

Other parameters:

-p/--pvalue: P-value based on a negative model in random introns, predefined thresholds in scores.py. 0.05 (default and recommended), 0.01, 0.005, none. If none is selected, all predicted ORFs will be printed.

-t/--threshold: Minimum predicted ORF length in amino acids. 24 (default) or 60.

-n/--orf_number: Select if the 'longest' (default) or 'all' ORFs per transcript are analysed.
```

**Example to run the test on cipher**:
```
python cipher.py -i test/sequences_human.fa -o test/test -s human -x tables/hsa_coding_to_intron_dicodon_usage.obj         
```

Online tool: http://evolutionarygenomics.imim.es/cipher



