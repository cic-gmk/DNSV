# coding=utf-8
# !/usr/bin/env python3
import os, re
import numpy as np
import pandas as pd
from SimpleCalculate import simpleStatistics

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
    #  if 'SVLEN' in sv_data['INFO'].iloc[0]: ?.*SVLEN=(?P<sv_len>-?[0-9]+)
    datGrab = re.compile("^.*SVTYPE=(?P<sv_type>[a-zA-Z]+).*$")
    data_info = datGrab.search(str(sv_data['INFO'].iloc[0])).groupdict()
    sv_type = data_info['sv_type']
    return sv_type

def svLen(sv_data):
    #  if 'SVLEN' in sv_data['INFO'].iloc[0]: ?.*SVLEN=(?P<sv_len>-?[0-9]+)
    datGrab = re.compile("^.*SVLEN=(?P<sv_len>-?[0-9]+).*$")
    if 'SVLEN' in str(sv_data['INFO'].iloc[0]):
        data_info = datGrab.search(str(sv_data['INFO'].iloc[0])).groupdict()
        sv_len = data_info['sv_len']
    else:
        # if the sv_type is not DEL or INS, we prefer to preserve it thus default sv_len 51 (>50).
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

def process_bar(i):
    num = i // 2
    if i == 100:
        process = "\r[%3s%%]: |%-50s|\n" % (i, '|' * num)
    else:
        process = "\r[%3s%%]: |%-50s|" % (i, '|' * num)
    print(process, end='', flush=True)
def judgeIfOverlap(start_1,end_1,start_2,end_2,sv_type,refdist,overlap_rate=0.5):
    #start_1 < end_1 && start_2 < end_2
    if sv_type == 'INS':
        return start_1 - refdist <= start_2 <= start_1 + refdist
    else:
        start_max = max(start_1,start_2)
        end_min = min(end_1,end_2)
        SV_1_range = end_1 - start_1
        SV_2_range = end_2 - start_2
        overlap_range = end_min - start_max
        return overlap_range/SV_1_range >=overlap_rate and overlap_range/SV_2_range >=overlap_rate
def getStartAndEnd(start,end):
    start,end = min(start,end),max(start,end)
    return start,end
def judgeNeighbour(bench_SVs,home_pos,compared_sv_chrom,compared_sv_start,compared_sv_end,compared_sv_type,typeignore,refdist,overlap_rate):
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
        the_TP = bench_SVs.xs(compared_sv_chrom).to_frame().T.iloc[[left_neighbour_origin_loc]]
    else:
        left_neighbour_end = svEnd(bench_SVs.xs(compared_sv_chrom).iloc[[left_neighbour_origin_loc]])
        left_neighbour_type = svType(bench_SVs.xs(compared_sv_chrom).iloc[[left_neighbour_origin_loc]])
        the_TP = bench_SVs.xs(compared_sv_chrom).iloc[[left_neighbour_origin_loc]]
    left_neighbour_start,left_neighbour_end = getStartAndEnd(pos_sort_list[left_neighbour_loc],left_neighbour_end)
    if judgeIfOverlap(compared_sv_start, compared_sv_end, left_neighbour_start,left_neighbour_end,compared_sv_type,refdist,overlap_rate):
        if typeignore == False:
            if left_neighbour_type == compared_sv_type:
                flag = 1
        else:
            flag = 1
    if flag == 1:
        return flag,the_TP
    elif flag == 0:
        if right_neighbour_loc is not None:
            right_neighbour_origin_loc = pos_list.index(pos_sort_list[right_neighbour_loc])
            if bench_SVs.xs(compared_sv_chrom)['POS'].shape == ():
                right_neighbour_end = svEnd(bench_SVs.xs(compared_sv_chrom).to_frame().T.iloc[[right_neighbour_origin_loc]])
                right_neighbour_type = svType(bench_SVs.xs(compared_sv_chrom).to_frame().T.iloc[[right_neighbour_origin_loc]])
                the_TP = bench_SVs.xs(compared_sv_chrom).to_frame().T.iloc[[right_neighbour_origin_loc]]
            else:
                right_neighbour_end = svEnd(bench_SVs.xs(compared_sv_chrom).iloc[[right_neighbour_origin_loc]])
                right_neighbour_type = svType(bench_SVs.xs(compared_sv_chrom).iloc[[right_neighbour_origin_loc]])
                the_TP = bench_SVs.xs(compared_sv_chrom).iloc[[right_neighbour_origin_loc]]
            right_neighbour_start,right_neighbour_end = getStartAndEnd(pos_sort_list[right_neighbour_loc],right_neighbour_end)
            if judgeIfOverlap(compared_sv_start, compared_sv_end, right_neighbour_start,right_neighbour_end,compared_sv_type,refdist,overlap_rate):
                if typeignore == False:
                    if right_neighbour_type == compared_sv_type:
                        flag = 1
                else:
                    flag = 1
    if flag == 1:
        return flag,the_TP
    else:
        return flag,None
def judgeIfSame(data_1,data_2,refdist,typeignore,overlap_rate,i):
    flag = 0
    # data_1_sv.shape[0]
    data_1_sv_pos = int(data_1['POS'][i])
    data_1_sv_chrom = data_1.index[i]
    data_1_sv_end = svEnd(data_1.iloc[[i]])
    data_1_sv_start,data_1_sv_end = getStartAndEnd(data_1_sv_pos, data_1_sv_end)
    data_1_sv_type = svType(data_1.iloc[[i]])
    # if the chrom is the same
    if data_1_sv_chrom in data_2.index:
        flag,the_TP = judgeNeighbour(data_2,data_1_sv_start,data_1_sv_chrom,data_1_sv_start,data_1_sv_end,data_1_sv_type,typeignore,refdist,overlap_rate)
        if data_1_sv_type != 'INS' and flag == 0:
            flag,the_TP = judgeNeighbour(data_2,data_1_sv_end,data_1_sv_chrom,data_1_sv_start,data_1_sv_end,data_1_sv_type,typeignore,refdist,overlap_rate)
   
    return flag,the_TP
def compareTP(file_name_1,file_name_2,out_dir,refdist,typeignore=False,overlap_rate=0.5):
    # dnsv
    if 'vcf' in file_name_1:
        data_1  = readvcf(file_name_1)
    else:
        data_1  = pd.read_csv(file_name_1,index_col='CHROM')
    # benchmark
    if 'vcf' in file_name_2:
        data_2  = readvcf(file_name_2)
    else:
        data_2  = pd.read_csv(file_name_2,index_col='CHROM')
    TP_sv = pd.DataFrame(columns=data_1.columns)
    TP_sv.index.name = 'CHROM'
    TP_benchmark_sv = pd.DataFrame(columns=data_2.columns)
    TP_benchmark_sv.index.name = 'CHROM'

    process_count = 0; process_path = data_1.shape[0]/100
    for i in range(data_1.shape[0]):
        #data_1.shape[0]
        if i >= process_path * process_count:
            process_bar(process_count+1)
            process_count = process_count + 1

        try:
#             data_1_sv_type = svType(data_1.iloc[[i]])
#             if data_1_sv_type not in ['INS','DEL']:
#                 continue
            flag,the_TP = judgeIfSame(data_1,data_2,refdist,typeignore,overlap_rate,i)
            if flag ==1:
                TP_sv = pd.concat([TP_sv, data_1.iloc[[i]]])
                TP_benchmark_sv = pd.concat([TP_benchmark_sv, the_TP])

        except:
            print(i)
            continue
    print(TP_sv)
    print(TP_benchmark_sv)
    
    statistics_output = simpleStatistics(TP_sv)
    TP_benchmark_statistics_output = simpleStatistics(TP_benchmark_sv)
    output_main = out_dir[:out_dir.rindex('.')]
    statistics_output_path = output_main+'_statistics.csv'
    TP_benchmark_statistics_output_path = output_main+'_benchmark_statistics.csv'
    TP_benchmark_output_path = output_main+'_benchmark.csv'
    
    TP_sv.to_csv(out_dir)
    TP_benchmark_sv.to_csv(TP_benchmark_output_path)
    statistics_output.to_csv(statistics_output_path)
    TP_benchmark_statistics_output.to_csv(TP_benchmark_statistics_output_path)
    
    print('Statistics Done!')

# if __name__ == '__main__':
    #AJson
    #pbsv
#     compareTP(file_name_1=r'pbsv_AJson_filtered_repeat_limit.csv',\
#            file_name_2=r'Tier1_pass_fixed.csv',\
#            out_dir = r'pbsv_AJson_filtered_repeat_limit_TP.csv',\
#            refdist=200,typeignore=False)
     