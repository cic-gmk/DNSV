# coding=utf-8
# !/usr/bin/env python3
import os, re
import numpy as np
import pandas as pd


def svLen(sv_data):
    data_grab = re.compile("^.*SVLEN=(?P<sv_len>-?[0-9]+).*$")
    if 'SVLEN' in str(sv_data['INFO'].iloc[0]):
        data_info = data_grab.search(sv_data['INFO'].iloc[0]).groupdict()
        sv_len = data_info['sv_len']
    else:
        # if the sv_type is not DEL, INS, DUP or INV, we prefer to preserve it thus default sv_len 51 (>50).
        sv_len = 51

    return int(sv_len)
def svType(sv_data):
    data_grab = re.compile("^.*SVTYPE=(?P<sv_type>[a-zA-Z]+).*$")
    if 'SVTYPE' in str(sv_data['INFO'].iloc[0]):
        data_info = data_grab.search(sv_data['INFO'].iloc[0]).groupdict()
        sv_type = data_info['sv_type']
    else:
        sv_type = 'None'
    return sv_type
def readvcf(file_name):
    count_num = 0
    with open(file_name,'r') as f1:
        for row in f1:
            if '#' in row:
                count_num = count_num + 1
#     print(count_num)
    rawData = pd.read_csv(file_name,skiprows=count_num-1,sep='\t')
    rawData = rawData.set_index('#CHROM')
    rawData.index.name = 'CHROM'
    # print(rawData.loc['chr1'])

    return rawData
def typeCalculate(file_name):
    if 'vcf' in file_name:
        sv_data = readvcf(file_name)
    else:
        sv_data = pd.read_csv(file_name)
#     print(sv_data)
#     dnsv_filter_data =pd.DataFrame(columns=dnsv_data.columns)
    sv_type_list = []

    for i in range(sv_data.shape[0]):
        print(i)
#         sv_len =svLen(sv_data.iloc[[i]])
#         if sv_len>10000:
        sv_type = svType(sv_data.iloc[[i]])
        sv_type_list.append(sv_type)
    sv_type_list = pd.Series(sv_type_list)
    print(sv_type_list.value_counts())
    return
def process_bar(i):
    num = i // 2
    if i == 100:
        process = "\r[%3s%%]: |%-50s|\n" % (i, '|' * num)
    else:
        process = "\r[%3s%%]: |%-50s|" % (i, '|' * num)
    print(process, end='', flush=True)
def calcultateImprecise(file_name):
    data = pd.read_csv(file_name)
    imprecise_ins = pd.DataFrame(columns=data.columns)
    imprecise_del = pd.DataFrame(columns=data.columns)
    imprecise = pd.DataFrame(columns=data.columns)

    process_count = 0; process_path = data.shape[0]/100
    for i in range(data.shape[0]):
        if i >= process_path * process_count:
            process_bar(process_count+1)
            process_count = process_count + 1
            
        sv_type =svType(data.iloc[[i]])
        if 'IMPRECISE' in data['INFO'].iloc[i]:
            imprecise = pd.concat([imprecise, data.iloc[[i]]])
            if sv_type == 'INS':
                imprecise_ins = pd.concat([imprecise_ins, data.iloc[[i]]])
            elif sv_type == 'DEL':
                imprecise_del  = pd.concat([imprecise_del , data.iloc[[i]]])

    print('ins',imprecise_ins)
    print('del',imprecise_del)
    print('all',imprecise)
#     deimprecise.to_csv(out_dir,index=None)
    return     

def filterImprecise(file_name,out_dir):
    data = pd.read_csv(file_name)
    deimprecise_ins = pd.DataFrame(columns=data.columns)
    deimprecise_del = pd.DataFrame(columns=data.columns)
    deimprecise = pd.DataFrame(columns=data.columns)

    process_count = 0; process_path = data.shape[0]/100
    for i in range(data.shape[0]):
        if i >= process_path * process_count:
            process_bar(process_count+1)
            process_count = process_count + 1
            
            sv_type =svType(data.iloc[[i]])
        if 'IMPRECISE' not in data['INFO'].iloc[i]:
            deimprecise = pd.concat([deimprecise, data.iloc[[i]]])
            if sv_type == 'INS':
                deimprecise_ins = pd.concat([deimprecise_ins, data.iloc[[i]]])
            elif sv_type == 'DEL':
                deimprecise_del  = pd.concat([deimprecise_del , data.iloc[[i]]])

    print('ins',deimprecise_ins)
    print('del',deimprecise_del)
    deimprecise.to_csv(out_dir,index=None)
    return 
def sizeChromStatistics(certain_type_data):
    # print(certain_type_data)

    # print(certain_type_data['CHROM'].value_counts())
    statistics_total = certain_type_data.shape[0]
    statistics_100bp = 0
    statistics_100bp_300bp = 0
    statistics_300bp_1kb = 0
    statistics_1kb = 0

    for i in range(certain_type_data.shape[0]):
        sv_len = abs(svLen(certain_type_data.iloc[[i]]))
        if sv_len < 100:
            statistics_100bp = statistics_100bp + 1
        elif 100<=sv_len<300:
            statistics_100bp_300bp = statistics_100bp_300bp + 1
        elif 300<=sv_len<1000:
            statistics_300bp_1kb = statistics_300bp_1kb + 1
        elif sv_len>=1000:
            statistics_1kb = statistics_1kb + 1
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
    statistics_list = [statistics_total,statistics_100bp,statistics_100bp_300bp,statistics_300bp_1kb,statistics_1kb]
    statistics_list.extend(chrom_part)

    return statistics_list
def simpleStatistics(file_name,out_dir=None):
    if 'vcf' in file_name:
        sv_data = readvcf(file_name)
    elif '.' in file_name:
        sv_data = pd.read_csv(file_name,index_col='CHROM')
    else:
        sv_data =file_name

    INS_data = pd.DataFrame(columns=sv_data.columns)
    DEL_data = pd.DataFrame(columns=sv_data.columns)
    other_type_data = pd.DataFrame(columns=sv_data.columns)

    statistics_output_index = ['INS', 'DEL', 'Other Types', 'Total']
    statistics_output_columns = ['Total', 'size<100bp', '100bp<=size<300bp', '300bp<=size<1kb','size>=1kb']
    statistics_output_columns.extend(['chr' + str(i) for i in range(1, 23)])
    statistics_output_columns.extend(['chrX','chrY','Other Chroms'])

    statistics_output = pd.DataFrame(index=statistics_output_index, columns=statistics_output_columns)

    for i in range(sv_data.shape[0]):
        # print(i)
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
    print(statistics_output)
    # statistics_output.to_csv(out_dir)
    return statistics_output
def lenStatistics(file_name):
    if 'vcf' in file_name:
        sv_data = readvcf(file_name)
    else:
        sv_data = pd.read_csv(file_name,index_col='CHROM')

    INS_data = pd.DataFrame(columns=sv_data.columns)
    DEL_data = pd.DataFrame(columns=sv_data.columns)
    other_type_data = pd.DataFrame(columns=sv_data.columns)

    count_over50bp = 0
    count_over30bp = 0
    count_under30bp = 0
    for i in range(sv_data.shape[0]):
        print(i)
        try:
            sv_len = svLen(sv_data.iloc[[i]])
            if sv_len>=50:
                count_over50bp = count_over50bp+1
            elif 30<=sv_len<50:
                count_over30bp = count_over30bp+1
            else:
                count_under30bp = count_under30bp+1
            
        except:
            print('Data Format Error Loc: %s ' % (i + 1))
            continue
    print('shape:',sv_data.shape[0])
    print('count_over50bp:',count_over50bp)
    print('count_over30bp:',count_over30bp)
    print('count_under30bp:',count_under30bp)
def svCoverage(sv_data):
    data_grab = re.compile("^.*:(?P<sv_DR>-?[0-9]+):(?P<sv_DV>-?[0-9]+)$")
#     print(sv_data.iloc[0,8])
    try:
        data_info = data_grab.search(sv_data.iloc[0,8]).groupdict()
        sv_DR = data_info['sv_DR']
        sv_DV = data_info['sv_DV']
    except:
        sv_DR = 0
        sv_DV = 0
#     print(sv_DR)
#     print(sv_DV)
    return sv_DR,sv_DV

def coverageDistribution(file_name,out_dir=None):
    if 'vcf' in file_name:
        sv_data = readvcf(file_name)
    else:
        sv_data = pd.read_csv(file_name,index_col='CHROM')
    DR_list = []
    DV_list = []
    DV_over25 = 0
    for i in range(sv_data.shape[0]):
#         print(sv_data.iloc[i][8])
        sv_DR,sv_DV = svCoverage(sv_data.iloc[[i]])
        DR_list.append(int(sv_DR))
        DV_list.append(int(sv_DV))
        if int(sv_DV) <= 20:
            DV_over25 += 1 
    DR_list = pd.DataFrame(DR_list)
    DV_list = pd.DataFrame(DV_list)
    print(DR_list.describe())
    print(DV_list.describe())
    print('DV_over25:',DV_over25)
    return
