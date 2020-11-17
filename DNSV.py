#coding=utf-8
#!/usr/bin/env python3
import os, re, argparse, sys
import numpy as np
import pandas as pd
import time

USAGE="""\
123
"""

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
        sv_type = 'NoType'

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

def judgeIfOverlap(start_1,end_1,start_2,end_2):
    #start_1 < end_1 && start_2 < end_2
    start_max = max(start_1,start_2)
    end_min = min(end_1,end_2)
    return start_max <= end_min
def getStartAndEnd(start,end):
    start,end = min(start,end),max(start,end)
    return start,end
def judgeNeighbour(bench_SVs,home_pos,compared_sv_chrom,compared_sv_start,compared_sv_end,compared_sv_type,typeignore):
    flag = 0
    if bench_SVs.xs(compared_sv_chrom)['POS'].shape == ():
        pos_list = list([bench_SVs.xs(compared_sv_chrom)['POS']])
    else:
        pos_list = list(bench_SVs.xs(compared_sv_chrom)['POS'])

    pos_list.append(home_pos)
    pos_sort_list = sorted(pos_list)
    home_loc = pos_sort_list.index(home_pos)
    left_neighbour_loc = home_loc - 1
    if home_loc + 1 < len(pos_sort_list):
        right_neighbour_loc = home_loc + 1
    else:
        right_neighbour_loc = None
    left_neighbour_origin_loc = pos_list.index(pos_sort_list[left_neighbour_loc])
    if bench_SVs.xs(compared_sv_chrom)['POS'].shape == ():
        left_neighbour_end = svEnd(bench_SVs.xs(compared_sv_chrom).to_frame().T.iloc[[left_neighbour_origin_loc]])
        left_neighbour_type = svType(bench_SVs.xs(compared_sv_chrom).to_frame().T.iloc[[left_neighbour_origin_loc]])
    else:
        left_neighbour_end = svEnd(bench_SVs.xs(compared_sv_chrom).iloc[[left_neighbour_origin_loc]])
        left_neighbour_type = svType(bench_SVs.xs(compared_sv_chrom).iloc[[left_neighbour_origin_loc]])
    left_neighbour_start,left_neighbour_end = getStartAndEnd(pos_sort_list[left_neighbour_loc],left_neighbour_end)
    if judgeIfOverlap(compared_sv_start, compared_sv_end, left_neighbour_start,left_neighbour_end):
        if typeignore == False:
            if left_neighbour_type == compared_sv_type:
                flag = 1
        else:
            flag = 1
    if flag == 0:
        if right_neighbour_loc is not None:
            right_neighbour_origin_loc = pos_list.index(pos_sort_list[right_neighbour_loc])
            if bench_SVs.xs(compared_sv_chrom)['POS'].shape == ():
                right_neighbour_end = svEnd(bench_SVs.xs(compared_sv_chrom).to_frame().T.iloc[[right_neighbour_origin_loc]])
                right_neighbour_type = svType(bench_SVs.xs(compared_sv_chrom).to_frame().T.iloc[[right_neighbour_origin_loc]])
            else:
                right_neighbour_end = svEnd(bench_SVs.xs(compared_sv_chrom).iloc[[right_neighbour_origin_loc]])
                right_neighbour_type = svType(bench_SVs.xs(compared_sv_chrom).iloc[[right_neighbour_origin_loc]])
            right_neighbour_start,right_neighbour_end = getStartAndEnd(pos_sort_list[right_neighbour_loc],right_neighbour_end)
            if judgeIfOverlap(compared_sv_start, compared_sv_end, right_neighbour_start,right_neighbour_end):
                if typeignore == False:
                    if right_neighbour_type == compared_sv_type:
                        flag = 1
                else:
                    flag = 1
    return flag

def judgeIfDenovo(father_SVs,mother_SVs,son_SVs,refdist,typeignore,i):
    flag = 0
    # son_sv.shape[0]
    son_sv_pos = int(son_SVs['POS'][i])
    son_sv_chrom = son_SVs.index[i]
    son_sv_end = svEnd(son_SVs.iloc[[i]])
    son_sv_start,son_sv_end = getStartAndEnd(son_sv_pos, son_sv_end)
    son_sv_type = svType(son_SVs.iloc[[i]])
    # if the chrom is the same
    if son_sv_chrom in father_SVs.index:
        #
        flag = judgeNeighbour(father_SVs,son_sv_start,son_sv_chrom,son_sv_start-refdist,son_sv_end+refdist,son_sv_type,typeignore)
        if flag == 0:
            flag = judgeNeighbour(father_SVs,son_sv_end,son_sv_chrom,son_sv_start-refdist,son_sv_end+refdist,son_sv_type,typeignore)
    if flag == 0:
        if son_sv_chrom in mother_SVs.index:
            flag = judgeNeighbour(mother_SVs,son_sv_start,son_sv_chrom,son_sv_start-refdist,son_sv_end+refdist,son_sv_type,typeignore)
            if flag == 0:
                flag = judgeNeighbour(mother_SVs,son_sv_end,son_sv_chrom,son_sv_start-refdist,son_sv_end+refdist,son_sv_type,typeignore)

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
def sizeChromStatistics(certain_type_data):
    # print(certain_type_data)

    # print(certain_type_data['CHROM'].value_counts())
    statistics_total = certain_type_data.shape[0]
    statistics_100bp = 0
    statistics_100bp_1kb = 0
    statistics_1kb_100kb = 0
    statistics_100kb = 0

    for i in range(certain_type_data.shape[0]):
        sv_len = abs(svLen(certain_type_data.iloc[[i]]))
        if sv_len < 100:
            statistics_100bp = statistics_100bp + 1
        elif 100<=sv_len<1000:
            statistics_100bp_1kb = statistics_100bp_1kb + 1
        elif 1000<=sv_len<1000000:
            statistics_1kb_100kb = statistics_1kb_100kb + 1
        elif sv_len>=1000000:
            statistics_100kb = statistics_100kb + 1
    #chr1 to chr22
    chrom_part = []
    for i in range(1,23):
        if 'chr'+str(i) in certain_type_data.index:
            chrom_part.append(certain_type_data.index.value_counts()['chr'+str(i)])
        else:
            chrom_part.append(0)
    #chrX chrY & Other Chroms
    if 'chrX' in certain_type_data.index:
        chrom_part.append(certain_type_data.index.value_counts()['chrX'])
    else:
        chrom_part.append(0)
    if 'chrY' in certain_type_data.index:
        chrom_part.append(certain_type_data.index.value_counts()['chrY'])
    else:
        chrom_part.append(0)
    chrom_part.append(statistics_total-sum(chrom_part))
    statistics_list = [statistics_total,statistics_100bp,statistics_100bp_1kb,statistics_1kb_100kb,statistics_100kb]
    statistics_list.extend(chrom_part)

    return statistics_list
def statistics(sv_data):
    # if 'vcf' in file_name:
    #     sv_data = readvcf(file_name)
    # else:
    #     sv_data = pd.read_csv(file_name,index_col='CHROM')

    INS_data = pd.DataFrame(columns=sv_data.columns)
    DEL_data = pd.DataFrame(columns=sv_data.columns)
    other_type_data = pd.DataFrame(columns=sv_data.columns)

    statistics_output_index = ['INS', 'DEL', 'Other Types', 'Total']
    statistics_output_columns = ['Total', 'size<100bp', '100bp<=size<1kb', '1kb<=size<100kb','size>=100kb']
    statistics_output_columns.extend(['chr' + str(i) for i in range(1, 23)])
    statistics_output_columns.extend(['chrX','chrY','Other Chroms'])

    statistics_output = pd.DataFrame(index=statistics_output_index, columns=statistics_output_columns)

    for i in range(sv_data.shape[0]):
        try:
            sv_type = svType(sv_data.iloc[[i]])
            if sv_type == 'INS':
                INS_data = pd.concat([INS_data, sv_data.iloc[[i]]])
            elif sv_type == 'DEL':
                DEL_data = pd.concat([DEL_data, sv_data.iloc[[i]]])
            else:
                other_type_data = pd.concat([other_type_data, sv_data.iloc[[i]]])
        except:
            print('Data Format Error Loc: %s ' % (i + 1))
            continue

    statistics_output.loc['INS']=sizeChromStatistics(INS_data)
    statistics_output.loc['DEL']=sizeChromStatistics(DEL_data)
    statistics_output.loc['Other Types'] = sizeChromStatistics(other_type_data)
    statistics_output.loc['Total'] = sizeChromStatistics(sv_data)

    return statistics_output

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
    thresg.add_argument("-r", "--refdist", type=int, default=500, \
                        help="Max reference location distance (%(default)dbp)")
    thresg.add_argument("-t", "--typeignore", type=bool, default=False, \
                        help="Variant types don't need to match to compare (False)")
    filteg = parser.add_argument_group("Filtering Arguments")
    filteg.add_argument("-p", "--precisionlimit", type=bool, default=True, \
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
    father_SVs = readvcf(args.father_SVs)
    mother_SVs = readvcf(args.mother_SVs)
    son_SVs = readvcf(args.son_SVs)

    dnsv = pd.DataFrame(columns=son_SVs.columns)
    dnsv.index.name = 'CHROM'

    print('DNSV Detection Start!')
    process_count = 0; process_path = son_SVs.shape[0]/100
    for i in range(son_SVs.shape[0]):
        if i >= process_path * process_count:
            process_bar(process_count+1)
            process_count = process_count + 1

        flag = judgeIfDenovo(father_SVs,mother_SVs,son_SVs,args.refdist,args.typeignore,i)
        if flag == 0:
            dnsv = pd.concat([dnsv, son_SVs.iloc[[i]]])
    print('DNSV Detection Done!')
    dnsv_filtered = dnsvFilter(dnsv=dnsv,precisionlimit=args.precisionlimit,sizemin=args.sizemin,typelimit=args.typelimit)
    dnsv_filtered.to_csv(args.output,index_label='CHROM')
    print('DNSV Filter Done!')
    if args.statistics is True:
        statistics_output = statistics(dnsv_filtered)
        output_main = args.output[:args.output.rindex('.')]
        statistics_output_path = output_main+'_statistics.csv'
        statistics_output.to_csv(statistics_output_path)
        print('DNSV Statistics Done!')

    return

if __name__ == '__main__':
    runmain(sys.argv[1:])
    # python DNSV.py F:/things/long_reads/SV工作资料/sniffles/sniffles_father.vcf  F:/things/long_reads/SV工作资料/sniffles/sniffles_mother.vcf F:/things/long_reads/SV工作资料/sniffles/sniffles_son.vcf -o son_dnsv.csv --statistics True
    # readvcf(r'F:\things\long_reads\SV工作资料\pbsv\pbsv_father.vcf')
#     data = pd.read_csv(r'C:\Users\gmkly\Desktop\SV_about\honey\pbhoneytailsresults\pbhoney_dnsv.csv',)
#     data = data.rename(columns={'index':'CHROM'})
#     print(data)
    