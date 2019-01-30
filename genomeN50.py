#!/usr/bin/env python

import sys,os
import argparse
from Bio import SeqIO

def store_seq(file):
    data = {}
    f = open(file)
    for seq_record in SeqIO.parse(f,"fasta"):
        data[seq_record.id] = seq_record.seq.tostring()
    f.close()
    return data

def split_scaffold(id_sca,seq_sca,mingap,data_ctg):
    first_seq,first_gap,end = 0,0,0
    check = False
    number = len(seq_sca)
    for i in xrange(number):
        item = seq_sca[i]
        if item == 'N' or item == 'n':
            if not check:
                first_gap = i
                check = True
        elif i == number - 1:
            id_ctg = id_sca + '_' + str(first_seq + 1 ) + '_' + str(i - first_seq + 1)
            data_ctg[id_ctg] = seq_sca[first_seq:i + 1]
        else:
            if check:
                if len(seq_sca[first_gap:i]) >= mingap:
                    #print mingap,len(seq_sca[first_gap:i])
                    end = i
                    id_ctg = id_sca + '_' + str(first_seq + 1) + '_' + str(first_gap - first_seq)
                    data_ctg[id_ctg] = seq_sca[first_seq:first_gap]
                    first_seq = i
                check = False

def extract_contig(data_sca,len_gap):
    data_ctg = {}
    for id_sca in data_sca:
        split_scaffold(id_sca,data_sca[id_sca],len_gap,data_ctg)
    return data_ctg

def total(data):
    data_length = []
    num_total,num_100,num_2k = 0,0,0
    seq_total,seq_100,seq_2k = '','',''
    for key in data:
        seq = data[key]
        length = len(seq)
        data_length.append(length)
        if length >= 100:
            num_100 += 1
            seq_100 += seq
            if length >= 2000:
                num_2k += 1
                seq_2k += seq
        num_total += 1
        seq_total += seq
    return data_length,[len(seq_total),len(seq_100),len(seq_2k)],[num_total,num_100,num_2k]

def sum_quality(total,data):
    length,longest,longest_num = 0,0,0
    data_length,data_number = [],[]
    data_check = {0:False,0.5:True,0.6:True,0.7:True,0.8:True,0.9:True}
    data_per = [0,0.5,0.6,0.7,0.8,0.9]
    number = len(data)
    data_sort = sorted(data,reverse=True)
    for i in xrange(number):
        if data_sort[i] >= longest:
            longest = data_sort[i]
            longest_num += 1
        length += data_sort[i]
        for j in xrange(1,6):
            per_now = data_per[j]
            per_fore = data_per[j -1]
            check_now = data_check[per_now]
            check_fore = data_check[per_fore]
            if length >= (total * per_now) and check_now and not check_fore:
                data_length.append(data_sort[i])
                data_number.append(i + 1)
                data_check[per_now] = False
    data_length.append(longest)
    data_number.append(longest_num)
    return data_length,data_number

def qualify(data):
    data_length,length_sum,number_sum = total(data)
    length_total = length_sum[0]
    sum_length,sum_number = sum_quality(length_total,data_length)
    result_length = sum_length + length_sum
    result_number = sum_number + number_sum
    return result_length,result_number

def result_out(file,length_sca,number_sca,length_ctg,number_ctg):
    out = '%13s\t%-15s\t%-15s\t%-13s\t%-13s\t\n' % ('Type','ScaffoldLength','ScaffoldNumber','ContigLength','ContigNumber')
    type_stat = ['N50','N60','N70','N80','N90','Longest','Total','Length>100bp','Length>2kb']
    for i in xrange(9):
        out += '%13s\t%-15d\t%-15d\t%13d\t%13d\t\n' % (type_stat[i],length_sca[i],number_sca[i],length_ctg[i],number_ctg[i])
    f = open(file,'w')
    f.write(out)
    f.close()

def main(file_sca,len_gap,file_result):
    data_sca = store_seq(file_sca)
    data_ctg = extract_contig(data_sca,len_gap)
    length_sca,number_sca = qualify(data_sca)
    #print 'ctg',data_ctg.keys()
    length_ctg,number_ctg = qualify(data_ctg)
    result_out(file_result,length_sca,number_sca,length_ctg,number_ctg)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Stat N50-N90 value of contigs and scaffolds.Editor:TY''')
    parser.add_argument('-s','--scaffold',metavar='fasta',required=True,help='Sequences of scaffold in *.fa/*.fasta format.')
    parser.add_argument('-m','--minlen',metavar='int',type=int,default=1,help='Define minimum length of gap.Default: 1')
    parser.add_argument('-r','--result',metavar='file',help='Define file storing results.Default:$genome_quality.xls')
    #parser.add_argument('-w','--workdir',metavar='directory',help='Define work directory.Default:$genome_work/')
    #parser.add_argument('-o','--outdir',metavar='directory',help='Define result directory.Default:$genome_out/')
    args = parser.parse_args()
    #print args
    dir_current = os.getcwd() + '/'
    name_genome = os.path.splitext(os.path.basename(args.scaffold))[0]
    if args.result:
        file_result = args.result
    else:
        file_result = dir_current + name_genome + '_result.xls'

    main(args.scaffold,args.minlen,file_result)
'''
    if not args.workdir:
        dir_work = dir_current + name_genome + '_work/'
    else:
        dir_work = args.workdir + '/'
    if not args.outdir:
        dir_result = dir_current + name_genome + '_out/'
    else:
        dir_result = args.outdir + '/'

    if os.path.exists(dir_work) or os.path.exists(dir_result):
        print >>sys.stderr,'Warnning:' + dir_work + ' or ' + dir_result + 'exists!\nPlease delete them or change your work directory and result directory.\n'
        exit()
    else:
        os.mkdir(dir_work)
        os.mkdir(dir_result)
'''
