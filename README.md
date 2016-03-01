# CIPHER
Online tool: http://evolutionarygenomics.imim.es/josepl/pol_web/index.php/primera/home

CIPHER is a program for the prediction of coding sequences, or open reading frames (ORFs), in transcripts. It is based on the differences in hexamer nucleotide frequencies in sets of known coding and non-coding sequences for different species (see Algorithm). It was used for the first time in the analysis of putative coding transcripts identified by ribosome profiling (Ruiz-Orera et al., 2014).  The user should upload a set of DNA sequences in FASTA format. Possible uses of the program are the prediction of putative proteins in a set of assembled transcripts or the detection of additional peptides produced from a given mRNA. The analysis is restricted to ORFs starting with an ATG codon (Methyonine) and ending with a STOP codon. Empirical p-value cutoffs have been precalculated for ORFs of different minimum length in various model species. 

MAIN EXECUTABLE: cipher.py  [options]

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i FASTA, --input=FASTA
                        Transcript file in FASTA format (required).
  -o OUT_FILE, --outfile=OUT_FILE
                        Output file (required). Tab separated text file: name
                        <tab> ORF number <tab> ORF position <tab ORF length
                        <tab> Transcript length <tab> Coding Score.
  -s SP, --specie=SP    Specie to be analyzed (required). Prebuilt thresholds
                        and codon matrices for: human, primate, mouse,
                        zebrafish, arabidopsis, yeast
  -x DICODONS, --table=DICODONS
                        Python object with dicodon scores (required)
  -p P_VALUE, --pvalue=P_VALUE
                        P-value: 0.05, 0.01, 0.005. default=0.05
  -t THRESHOLD, --threshold=THRESHOLD
                        Minimum ORF length threshold in amino acids.
                        default=24
  -n ORF_S, --orf_number=ORF_S
                        Predict longest or all ORFs. default=longest

                        
FILE TO CREATE DICODONS TABLES: dicodon_matrices.py  [options]

Options:
  -h, --help            show this help message and exit
  -c C, --coding=C      Coding file: ORFs in nucleotide FASTA format
                        (required).
  -n NC, --noncoding=NC
                        Non-coding file: ORFs in nucleotide FASTA format
                        (required).
  -o OUT_FILE, --outfile=OUT_FILE
                        Output file, written in 'tables' directory (required).


scores.py contains a function with computed thresholds for different ORFs and p-values.


Example to run the test on cipher:

python cipher.py -i tutorial/sequences_human.fa -o test/test -s human -x tables/hsa_coding_to_intron_dicodon_usage.obj

References:
Long non-coding RNAs as a source of new peptides Ruiz-Orera, J., Messeguer, X., Subirana, J. A., & Alba, M. M. (2014). Long non-coding RNAs as a source of new peptides. eLife, 3, 1 ‒ 24.



