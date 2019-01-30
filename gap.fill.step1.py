#!/usr/bin/python -w
import sys
sys.path.insert(0,"/home/huj/.pip")
import re
import pysam
from itertools import product
from Bio.Seq import Seq

def cliplen(cigartuples):
    clip5 = cigartuples[0][1] if cigartuples[0][0] in [4, 5] else 0
    clip3 = cigartuples[-1][1] if cigartuples[-1][0] in [4, 5] else 0
    return clip5 + clip3

def determine_align_position(i):
	# g = re.findall(r'(\d+)(\D)', cigar)
	# target = ref ; query = reads 
	lable = True
	tstart = tend = qstart = qend = 0
	M_number = len(re.findall(r'M', i.cigarstring))
	for op, position in i.cigar:
		if lable:	
			if op == 3:
				tstart += position
			elif op == 0:
				lable = False
				tend += tstart + position
				qend += qstart + position
				M_number -= 1
			elif op == 1 or op == 4 or op == 5 or op == 6:
				qstart += position
			elif op == 2:
				tstart += 1
			else:
				print >>sys.stderr, i.cigarstring, 0, op
		elif M_number:
			if op == 1:
				qend += position
			elif op == 2:
				tend += position
			elif op == 0:
				qend += position
				tend += position
				M_number -= 1
			else:
				print >>sys.stderr, i.cigarstring, 1, op
	return (tstart, tend, qstart, qend)

def gap_positions(query_name, qstart, qend, number, strand,):
	gap_info = query_name.split("_")
	gap_chr = gap_info[0]
	p_start, p_end, m_start, m_end = map(int, gap_info[1:])
	if number == 0:
		if strand == 0:
			qstart = p_start + qstart
			qend = p_start + qend -1
		else:
			temp = p_end - qend + 1
			qend = p_end - qstart
			qstart = temp

	else:
		if strand == 0:
			qstart = m_start + qstart
			qend = m_start + qend - 1
		else:
			temp = m_end - qend + 1
			qend = m_end - qstart
			qstart = temp
	return qstart, qend



def read_sam_file(in_file):
	result = {}
	count = {}
	for i in pysam.AlignmentFile(in_file, "r"):
		read_mapped = 0 if i.flag & 4 else 1 #0:un 1:map
		mate_mapped = 0 if i.flag & 8 else 1 #0:un 1:map
		read_number = 0 if i.flag & 64 else 1 # 0:first reads 1:second reads
 		read_strand = 1 if i.flag & 16 else 0 # 0: + 1: -
 		if i.query_name not in result:
			result[i.query_name] = {}
			count[i.query_name] = 0
		if not read_mapped:
			continue
		if cliplen(i.cigartuples)/float(i.infer_read_length()) > 0.8:
			continue
		if read_number not in result[i.query_name]:
			result[i.query_name][read_number] = {}
		if i.reference_name not in result[i.query_name][read_number]:
			result[i.query_name][read_number][i.reference_name] = []
		(tstart, tend, qstart, qend) = determine_align_position(i)
		tstart += i.reference_start + 1
		tend += tstart - 1
		# print tstart, tend, qstart, qend, i.reference_start, i.reference_end, i.query_alignment_start, i.query_alignment_end, i.cigarstring
		qstart, qend = gap_positions(i.query_name, qstart, qend, read_number, read_strand) 
		result[i.query_name][read_number][i.reference_name].append([qstart, qend, tstart, tend,read_strand]) # target: ref, query: reads
		count[i.query_name] += 1
	return result, count

def check_gap_fill(k):
	for query_name in k:
		closeList = []
		ErrorOrient = False
		gap_info = query_name.split("_")
		gap_chr = gap_info[0]
		p_start, p_end, m_start, m_end = map(int, gap_info[1:])
		if len(k[query_name]) == 0:
			print "\t".join(map(str, ("unmap", gap_chr, p_end, m_start, m_start - p_end + 1 )))
		elif len(k[query_name]) == 1:
			print "\t".join(map(str, ("oneside", gap_chr, p_end, m_start, m_start - p_end + 1 )))  
		else:
			p_chr = k[query_name][0].keys()
			m_chr = k[query_name][1].keys()
			chrs = set(p_chr + m_chr)
			for i in chrs:
				if i in p_chr and i in m_chr:
					for m, n in product(k[query_name][0][i], k[query_name][1][i]): # qstart, qend, tstart, tend,read_strand
							if (m[-1] == n[-1] == 0 and m[3] <= n[2]) or (m[-1] == n[-1] == 1 and n[3] <= m[2]):
								closeList.append(m + n + [i])
							else:
								ErrorOrient = True
			if closeList:
				# for i in closeList:
				# 	print i, count[query_name], query_name
				out = sorted(closeList, key = lambda x:abs(abs(x[7] - x[3]) - abs(x[5] - x[1])))[0] #select the nearest gap size
				if abs(abs(out[7] - out[3]) - abs(out[5] - out[1])) > 500000 and (count[query_name] >= 5 \
					or (len(closeList)) != 1):continue # filter the super difference gap size between except and true
				if out[4] == 0:
				#type chr pstart pend plen truestart trueend truelen refchr refstart refend strand
					print "\t".join(map(str, ("Closed", gap_chr, p_end, m_start, m_start - p_end + 1, out[1], out[5], out[5] - out[1] + 1, out[-1], out[3], out[7], out[4] )))
				else:
					print "\t".join(map(str, ("Closed", gap_chr, p_end, m_start, m_start - p_end + 1, out[1], out[5], out[5] - out[1] + 1, out[-1], out[8], out[2], out[4] )))
			elif ErrorOrient:
				print "\t".join(map(str, ("ErrorOrient", gap_chr, p_end, m_start, m_start - p_end + 1 )))
			else:
				print "\t".join(map(str, ("Trans", gap_chr, p_end, m_start, m_start - p_end + 1 )))






alignSamFile = sys.argv[1]
align_reads, count = read_sam_file(alignSamFile)
check_gap_fill(align_reads)

