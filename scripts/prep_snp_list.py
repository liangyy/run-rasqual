import re
import pandas as pd

def parse_varid(varid, complete=False):
    o1 = []
    o2 = []
    o3 = []
    if complete is True:
        o4 = []
        o5 = []
    for v in varid:
        tmp = v.split('_')
        cc = re.sub('chr', '', tmp[0])
        pos = tmp[1]
        o1.append(cc)
        o2.append(pos)
        o3.append(v)
        if complete is True:
            o4.append(tmp[2])
            o5.append(tmp[3])
    if complete is False:
        return pd.DataFrame({'chr': o1, 'pos': o2, 'snp_id': o3})
    else:
        return pd.DataFrame({'#CHROM': o1, 'POS': o2, 'ID': o3, 'REF': o4, 'ALT': o5})

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='prep_snp_list.py', description='''
        Prepare SNP list.
    ''')
    parser.add_argument('--genotype', help='''
        genotype in parquet
    ''')
    parser.add_argument('--output', help='''
        output TSV.GZ
    ''')
    args = parser.parse_args()
 
    import logging, time, sys, os
    # configing util
    logging.basicConfig(
        level = logging.INFO, 
        stream = sys.stderr, 
        format = '%(asctime)s  %(message)s',
        datefmt = '%Y-%m-%d %I:%M:%S %p'
    )
    
    
    df = pd.read_parquet(args.genotype, columns='index')
    df = df.index.to_list()
    parse_varid(df).to_csv(args.output, 
    sep='\t', header=True, index=False, compression='gzip')
    
