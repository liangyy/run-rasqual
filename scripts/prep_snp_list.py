import re
import pandas as pd

def parse_varid(varid):
    o1 = []
    o2 = []
    o3 = []
    for v in varid:
        tmp = v.split('_')
        cc = re.sub('chr', '', tmp[0])
        pos = tmp[1]
        o1.append(cc)
        o2.append(pos)
        o3.append(v)
    return pd.DataFrame({'chr': cc, 'pos': pos, 'snp_id': v})

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
    