# Genome-data-of-Hanfu-apple
 FASTA files of chromosomes and genes, gff files for gene models, scripts.

## Genome data
1. HFTH1.genome.fa
The assembled genome sequences (Pseudo-chromosomes) with fasta format of 'Hanfu' apple.

2. HFTH1.gene.gff3.gz
Gff files for gene models

3. HFTH1.gene.cds.fa.gz
Coding sequences for gene models

4. HFTH1.gene.pep.fa.gz
Protein sequences for gene models.

## Script:
1. stat N50-N90 value of contigs and scaffolds: genomeN50.py

2. repeat annotation with Repeatmodeler: Run_DenovoRepeatFinder.pl 

3. repeat annotation with RepeatMasker, RepeatProteinMask, Trf: Run_Repeat_Protein_Mask_Trf.pl

4. mask repeat region with N: Run_remask.pl  

5. stat Repeat results: Run_RepeatStatistic.py

6. gap filling: gap.fill.sh gap.fill.step1.py gap.fill.step2.py