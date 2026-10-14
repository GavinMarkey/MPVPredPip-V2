(Written by setup/init_input.py from config.yaml. It is rewritten if the
 configuration changes, so put your own notes somewhere else.)

INPUT FILES
===========

Put your protein sequences directly in this folder:  0 - Input/
Subfolders are ignored.  Accepted extensions:  .fasta  .fa  .faa

1. WHAT TO PUT HERE
   One protein sequence per file, downloaded from NCBI Protein in FASTA format.
   Amino acids only - nucleotide FASTA is rejected. Minimum 50 residues.

2. THE HEADER LINE IS WHAT THE PIPELINE READS
       >ACCESSION protein name [strain]
   for example
       >NP_066246.1 spike glycoprotein [Zaire ebolavirus]

   accession  =  the first word
   protein    =  the text before the [ ... ]
   strain     =  the text inside the [ ... ]

   NCBI normally puts the SPECIES (its organism name) in the brackets, so by
   default each species counts as one "strain". Files downloaded from NCBI
   already look like this, so you rarely need to edit them.

3. NAME THE FILE
       <ACCESSION> <protein name> - <strain>.fasta
   for example
       NP_066246.1 spike glycoprotein - Zaire ebolavirus.fasta

   The filename is only a FALLBACK. It is used when a header has no [strain] or
   no protein name, and it needs a space on BOTH sides of the dash. When the
   header is complete the filename can be anything - but a consistent name keeps
   this folder readable and keeps the fallback working.

4. THINGS THAT STOP THE RUN AT STAGE 1
   * Every (protein, strain) pair must be unique. Two Zaire spike sequences
     with identical [bracket text] are rejected as duplicates. To supply
     several isolates of one species, make the bracketed text differ:
         [Zaire ebolavirus Makona]      [Zaire ebolavirus Mayinga]
   * A sequence that is nucleotide, empty, or shorter than 50 residues.

5. HOW A FILE IS ASSIGNED TO A PROTEIN
   The protein name is matched against sequences.protein_aliases in config.yaml.
   Currently configured:
       "spike glycoprotein"  ->  Spike
       "glycoprotein"  ->  Spike
       "nucleoprotein"  ->  Nucleoprotein

   Matching is by substring, longest alias first, so a name that merely
   CONTAINS an alias gets that label. The Ebola GP gene also encodes a "secreted
   glycoprotein"; that name contains "glycoprotein" and would be labelled
   Spike. A name matching no alias becomes a protein of its own (CamelCased),
   so add an alias for every protein you analyse.

CHECK YOUR FILES BEFORE RUNNING
   python setup/init_input.py

THEN RUN
   snakemake --use-conda --cores 8
