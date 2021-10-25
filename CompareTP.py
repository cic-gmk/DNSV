# coding=utf-8
# !/usr/bin/env python3
import os, re
import numpy as np
import pandas as pd
from ReadUtils import readFile,svType,svLen,svEnd,processBar
from SimpleCalculate import simpleStatistics


def judgeIfOverlap(start_1,end_1,start_2,end_2,sv_type,refdist,overlap_rate=0.5):
    #start_1 < end_1 && start_2 < end_2
    
    if sv_type == 'INS':
        return start_1 - refdist <= start_2 <= start_1 + refdist or \
                start_1 - refdist <= end_2 <= start_1 + refdist
    else:
        start_max = max(start_1,start_2)
        end_min = min(end_1,end_2)
        overlap_range = end_min - start_max
        SV_1_range = end_1 - start_1
        SV_2_range = end_2 - start_2
        return overlap_range/SV_1_range >=overlap_rate and overlap_range/SV_2_range >=overlap_rate

def getStartAndEnd(start,end):
    start,end = min(start,end),max(start,end)
    return start,end

def binarySearch(bench_df,home_pos,low,high):
    
    left = int(low + (high - low)/2)
    right = left + 1
    if left == 0 or right == bench_df.shape[0] - 1:
        return left,right
    left_pos = bench_df['POS'].iloc[left]
    right_pos = bench_df['POS'].iloc[right]
    if left_pos <= home_pos <= right_pos:
        return left,right
    elif home_pos < left_pos:
        return binarySearch(bench_df, home_pos, low, left)
    else:
        return binarySearch(bench_df, home_pos, right, high)
        
def judgeNeighbour(bench_df,home_pos,compared_sv_chrom,compared_sv_start,compared_sv_end,compared_sv_type,typeignore,refdist,overlap_rate):
    flag = 0
    
    if (bench_df.shape == ()):
        if (home_pos < bench_df['POS'].iloc[0]):
            left_neighbour_loc,right_neighbour_loc = None,0
        else:
            left_neighbour_loc,right_neighbour_loc = 0,None
    else:
        left_neighbour_loc,right_neighbour_loc = binarySearch(bench_df, home_pos, 0, bench_df.shape[0]-1)
    if left_neighbour_loc is not None:
        if bench_df['POS'].shape == ():
            left_neighbour_end = svEnd(bench_df.to_frame().T.iloc[[left_neighbour_loc]])
            left_neighbour_type = svType(bench_df.to_frame().T.iloc[[left_neighbour_loc]])           
            the_TP = bench_df.to_frame().T.iloc[[left_neighbour_loc]]
        else:
            left_neighbour_end = svEnd(bench_df.iloc[[left_neighbour_loc]])
            left_neighbour_type = svType(bench_df.iloc[[left_neighbour_loc]])
            the_TP = bench_df.iloc[[left_neighbour_loc]]
        left_neighbour_start,left_neighbour_end = getStartAndEnd(bench_df['POS'].iloc[left_neighbour_loc],left_neighbour_end)
        if judgeIfOverlap(compared_sv_start, compared_sv_end, left_neighbour_start,left_neighbour_end,compared_sv_type,refdist,overlap_rate):
            if typeignore == False:
                if left_neighbour_type == compared_sv_type:
                    flag = 1
            else:
                flag = 1
        if flag == 1:
            return flag,the_TP
        
    if right_neighbour_loc is not None:
        if bench_df['POS'].shape == ():
            right_neighbour_end = svEnd(bench_df.to_frame().T.iloc[[right_neighbour_loc]])
            right_neighbour_type = svType(bench_df.to_frame().T.iloc[[right_neighbour_loc]])
            the_TP = bench_df.to_frame().T.iloc[[right_neighbour_loc]]
        else:
            right_neighbour_end = svEnd(bench_df.iloc[[right_neighbour_loc]])
            right_neighbour_type = svType(bench_df.iloc[[right_neighbour_loc]])
            the_TP = bench_df.iloc[[right_neighbour_loc]]
        right_neighbour_start,right_neighbour_end = getStartAndEnd(bench_df['POS'].iloc[right_neighbour_loc],right_neighbour_end)
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
    data_1_sv_chrom = str(data_1.index[i])
    data_1_sv_end = svEnd(data_1.iloc[[i]])
    data_1_sv_start,data_1_sv_end = getStartAndEnd(data_1_sv_pos, data_1_sv_end)
    data_1_sv_type = svType(data_1.iloc[[i]])
    # if the chrom is the same
    
    
    if data_1_sv_chrom in data_2.index:
        global bench_dict
        bench_df = bench_dict[data_1_sv_chrom]
        flag,the_TP = judgeNeighbour(bench_df,data_1_sv_start,data_1_sv_chrom,data_1_sv_start,data_1_sv_end,data_1_sv_type,typeignore,refdist,overlap_rate)
        if flag == 0 and (data_1_sv_type not in ['INS','None']):
            flag,the_TP = judgeNeighbour(bench_df,data_1_sv_end,data_1_sv_chrom,data_1_sv_start,data_1_sv_end,data_1_sv_type,typeignore,refdist,overlap_rate)
    return flag,the_TP
    

def preProcessData(bench_data):
    global bench_dict
    bench_dict = {}
    for chrom in bench_data.index:
        bench_dict[chrom] = bench_data.xs(chrom).sort_values(by = 'POS')    

def compareTP(file_name_1,file_name_2,out_dir,refdist,typeignore=False,overlap_rate=0.5):
    # SV_1
    data_1  = readFile(file_name_1)
    # benchmark
    data_2  = readFile(file_name_2)
    
    print("Initialization Start!")
    preProcessData(data_2)
    print("Initialization Done!")
    
    TP_sv = pd.DataFrame(columns=data_1.columns)
    TP_sv.index.name = 'CHROM'
    TP_benchmark_sv = pd.DataFrame(columns=data_2.columns)
    TP_benchmark_sv.index.name = 'CHROM'

    process_count = 0; process_path = data_1.shape[0]/100
    for i in range(data_1.shape[0]):
        #data_1.shape[0]
        if i >= process_path * process_count:
            processBar(process_count)
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
    processBar(100)
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
     
