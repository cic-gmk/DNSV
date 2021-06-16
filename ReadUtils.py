
import pandas as pd
import numpy as np
import re
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
def processBar(i):
    num = i // 2
    if i == 100:
        process = "\r[%3s%%]: |%-50s|\n" % (i, '|' * num)
    else:
        process = "\r[%3s%%]: |%-50s|" % (i, '|' * num)
    print(process, end='', flush=True)