#coding=utf-8
#!/usr/bin/env python3
import os, re, argparse, sys
import numpy as np
import pandas as pd
from SimpleCalculate import simpleStatistics
from CompareOverlap import judgeNeighbour,getStartAndEnd

USAGE = 'Usage: python DNSV.py [options] father.vcf mother.vcf son.vcf -o dnsv.csv'
def readFile(file_name):
    if 'vcf' in file_name:
        data = readvcf(file_name)
    else:
        data = pd.read_csv(file_name,index_col='CHROM')
    return data

def readvcf(file_name):
    count_num = 0
    with open(file_name,'r') as f1:
        for row in f1:
            if '#' in row:
                count_num = count_num + 1
            else:
                break
#     print(count_num)
    raw_data = pd.read_csv(file_name,skiprows=count_num-1,sep='\t')
    raw_data = raw_data.set_index('#CHROM')
    raw_data.index.name = 'CHROM'

    return raw_data


def svType(sv_data):
    data_grab = re.compile("^.*SVTYPE=(?P<sv_type>[a-zA-Z]+).*$")
    if 'SVTYPE' in sv_data['INFO'].iloc[0]:
        data_info = data_grab.search(sv_data['INFO'].iloc[0]).groupdict()
        sv_type = data_info['sv_type']
    else:
        sv_type = 'None'

    return sv_type

def svLen(sv_data):
    data_grab = re.compile("^.*SVLEN=(?P<sv_len>-?[0-9]+).*$")
    if 'SVLEN' in sv_data['INFO'].iloc[0]:
        data_info = data_grab.search(sv_data['INFO'].iloc[0]).groupdict()
        sv_len = data_info['sv_len']
    else:
        # if the sv_type is not DEL, INS, DUP or INV, we prefer to preserve it thus default sv_len 51 (>50).
        sv_len = 51

    return int(sv_len)
def svEnd(sv_data):
    data_grab = re.compile("^.*END=(?P<sv_end>-?[0-9]+).*$")
    if 'END' in sv_data['INFO'].iloc[0]:
        data_info = data_grab.search(sv_data['INFO'].iloc[0]).groupdict()
        sv_end = data_info['sv_end']
    else:
        sv_end = sv_data['POS']
    return int(sv_end)


def judgeIfDenovo(father_SVs,mother_SVs,son_SVs,refdist,typeignore,overlap_rate,i):
    flag = 0
    # son_sv.shape[0]
    son_sv_pos = int(son_SVs['POS'][i])
    son_sv_chrom = son_SVs.index[i]
    son_sv_end = svEnd(son_SVs.iloc[[i]])
    son_sv_start,son_sv_end = getStartAndEnd(son_sv_pos, son_sv_end)
    son_sv_type = svType(son_SVs.iloc[[i]])
    # if the chrom is the same
    if son_sv_chrom in father_SVs.index:
        flag = judgeNeighbour(father_SVs,son_sv_start,son_sv_chrom,son_sv_start,son_sv_end,son_sv_type,typeignore,refdist,overlap_rate)
        if son_sv_type != 'INS' and flag == 0:
            flag = judgeNeighbour(father_SVs,son_sv_end,son_sv_chrom,son_sv_start,son_sv_end,son_sv_type,typeignore,refdist,overlap_rate)
    if flag == 0:
        if son_sv_chrom in mother_SVs.index:
            flag = judgeNeighbour(mother_SVs,son_sv_start,son_sv_chrom,son_sv_start,son_sv_end,son_sv_type,typeignore,refdist,overlap_rate)
            if son_sv_type != 'INS' and flag == 0:
                flag = judgeNeighbour(mother_SVs,son_sv_end,son_sv_chrom,son_sv_start,son_sv_end,son_sv_type,typeignore,refdist,overlap_rate)

    return flag

def dnsvFilter(dnsv,precisionlimit=True,sizemin=50,typelimit=True):
    # print(dnsv_data)
    dnsv_filtered =pd.DataFrame(columns=dnsv.columns)

    for i in range(dnsv.shape[0]):
        try:
            #dnsv_data.shape[0]
            # print(dnsv_data.iloc[[i]]['INFO'].iloc[0])
            if precisionlimit:
                if 'IMPRECISE' in dnsv['INFO'].iloc[i]:
                    continue
            sv_type = svType(dnsv.iloc[[i]])
            if typelimit:
                if  sv_type not in ['INS','DEL','DUP','INV']:
                    continue
            if sizemin is not None:
                if sv_type in ['INS','DEL','DUP','INV']:
                    if abs(svLen(dnsv.iloc[[i]])) < sizemin:
                        continue
            dnsv_filtered =pd.concat([dnsv_filtered, dnsv.iloc[[i]]])
        except:
            print('Data Format Error Loc: %s ' % (i+1))
            continue
    # dnsv_filtered_data.to_csv(out_dir)
    return dnsv_filtered


def process_bar(i):
    num = i // 2
    if i == 100:
        process = "\r[%3s%%]: |%-50s|\n" % (i, '|' * num)
    else:
        process = "\r[%3s%%]: |%-50s|" % (i, '|' * num)
    print(process, end='', flush=True)

def parseArgs(argv):
    parser = argparse.ArgumentParser(prog="DNSV.py", description=USAGE, \
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("father_SVs", type=str, \
                        help="Input father's SVs")
    parser.add_argument("mother_SVs", type=str, \
                        help="Input mother's SVs")
    parser.add_argument("son_SVs", type=str, \
                        help="Input son's SVs")
    parser.add_argument("-o", "--output", type=str, default=None, \
                        help="Output Name (son_sv_DNSV.csv)")
    parser.add_argument("--statistics", type=bool, default=False, \
                        help="Whether the statistics of DNSVs is needed (False)")
    thresg = parser.add_argument_group("Comparison Threshold Arguments")
    thresg.add_argument("-r", "--refdist", type=int, default=200, \
                        help="Max reference location distance (%(default)dbp)")
    thresg.add_argument("-t", "--typeignore", type=bool, default=False, \
                        help="Variant types don't need to match to compare (False)")
    thresg.add_argument("-O", "--overlap_rate", type=float, default=0.5, \
                        help="Reciprocal overlaps with the reference SVs (0.5)")
    filteg = parser.add_argument_group("Filtering Arguments")
    filteg.add_argument("-p", "--precisionlimit", type=bool, default=False, \
                        help="Limit the precision description in INFO column")
    filteg.add_argument("-s", "--sizemin", type=int, default=50, \
                        help="Minimum DNSV size (%(default)dbp)")
    filteg.add_argument("-l", "--typelimit", type=bool, default=True, \
                        help="Limit DNSV types in ['INS','DEL','DUP','INV']")

    args = parser.parse_args(argv)

    #     setupLogging(args.debug)

    if args.output is None:
        if '.' in args.son_sv:
            output_main = args.son_sv[:args.son_sv.rindex('.')]
        else:
            output_main = args.son_sv
        args.output = output_main + "_DNSV.csv"

    return args

def runmain(argv):
    # father_sv, vcf_mother, vcf_son, out_dir, refdist, typeignore = False
    args = parseArgs(argv)

    father_SVs  = readFile(args.father_SVs)
    mother_SVs  = readFile(args.mother_SVs)
    son_SVs  = readFile(args.son_SVs)

    dnsv = pd.DataFrame(columns=son_SVs.columns)
    dnsv.index.name = 'CHROM'

    print('DNSV Detection Start!')
    process_count = 0; process_path = son_SVs.shape[0]/100
    for i in range(son_SVs.shape[0]):

        if i >= process_path * process_count:
            process_bar(process_count+1)
            process_count = process_count + 1
        try:
            flag = judgeIfDenovo(father_SVs,mother_SVs,son_SVs,args.refdist,args.typeignore,args.overlap_rate,i)
            if flag == 0:
                dnsv = pd.concat([dnsv, son_SVs.iloc[[i]]])
        except:
            print(i)
            continue
    print('DNSV Detection Done!')
    dnsv_filtered = dnsvFilter(dnsv=dnsv,precisionlimit=args.precisionlimit,sizemin=args.sizemin,typelimit=args.typelimit)
    dnsv_filtered.to_csv(args.output,index_label='CHROM')
    print('DNSV Filter Done!')
    if args.statistics is True:
        statistics_output = simpleStatistics(dnsv_filtered)
        output_main = args.output[:args.output.rindex('.')]
        statistics_output_path = output_main+'_statistics.csv'
        statistics_output.to_csv(statistics_output_path)
        print('DNSV Statistics Done!')

    return

if __name__ == '__main__':
    runmain(sys.argv[1:])
#python DNSV.py F:\things\long_reads\SV工作资料\AJtrio\sniffles\docs\father_test.csv F:\things\long_reads\SV工作资料\AJtrio\sniffles\docs\mother_test.csv F:\things\long_reads\SV工作资料\AJtrio\sniffles\docs\son_test.csv -o F:\things\long_reads\SV工作资料\AJtrio\sniffles\docs\dnsv_test.csv --statistics True