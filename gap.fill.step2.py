#!/usr/bin/env python

import sys
import re
from filehandler import Fasta
from Bio.Seq import Seq


def trans(seqs):
    return Seq(seqs).reverse_complement()

def read_gap_file(infile):
	result = {}
	with open(infile) as IN:
		for line in IN:
			lines = line.strip().split()
			if lines[0] != "Closed":
				continue
			if lines[1] not in result:
				result[lines[1]] = []
			result[lines[1]].append([int(i) for i in lines[2:] if i.isdigit()] + [lines[8]])
	return result



#python fill.gap.py *.gap fill.fasta gap.fasta
gapIn = read_gap_file(sys.argv[1])
target_seq = Fasta(sys.argv[2]).readRecords()


for sca, seqs in Fasta(sys.argv[3]).readRecord():
	if sca in gapIn:
		seq_new = ''
		preend = 0
		gap_sort = sorted(gapIn[sca], key=lambda x:x[0])
		for k in gap_sort:
			# Closed  Chr08   7869989 7906716 36728   7869989 7906716 36728   Scaffold25_1_18010644   7398278 7414583 0
			gap_start, gap_end, fill_start, fill_end, strand, fill_chr = k[3:5] + k[6:]
			if strand == "0":

				gap = seqs[gap_start -5:gap_end + 5]
				# target = target_seq[fill_chr][fill_start-5:fill_end + 5]
			else:
				gap = seqs[gap_start -5:gap_end + 5]
			# 	target = trans(target_seq[fill_chr][fill_start-6:fill_end + 4])
			# print ">" + sca + "_" + str(gap_start -5) +str("_") + str(gap_end + 5)
			# print target
			seq_new += seqs[preend:gap_start -5] + target_seq[sca + "_" + str(gap_start -5) +str("_") + str(gap_end + 5) + "_pilon"]
			preend  =  gap_end + 5
		seq_new += seqs[preend:]
		seqs = seq_new
	print ">" + sca + "\n" + seqs