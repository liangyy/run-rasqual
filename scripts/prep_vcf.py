import pandas as pd
from tqdm import tqdm
from prep_snp_list import parse_varid
 
def sample2indiv(str_):
    o = []
    for i in str_:
        tmp = str_.split('-')
        o.append('-'.join(tmp[:2]))
    return o

def construct_vcf(ref, alt, h1, h2):
    df_template = pd.DataFrame({'VARIANT_ID': h1.index.to_list()})
    ref = pd.merge(df_template, ref, on='VARIANT_ID', how='left').fillna(0)
    alt = pd.merge(df_template, alt, on='VARIANT_ID', how='left').fillna(0)
    mat = _get_gt_and_as(
        ref.drop(columns='VARIANT_ID').values, 
        alt.drop(columns='VARIANT_ID').values, 
        h1.values, h2.values
    )
    mat = pd.DataFrame(mat)
    mat.columns = h1.columns
    df_vcf = parse_varid(ref.VARIANT_ID.to_list())
    df_vcf['QUAL'] = 100
    df_vcf['FILTER'] = 'PASS'
    df_vcf['INFO'] = 'RSQ=0.99'
    df_vcf['FORMAT'] = 'GT:AS'
    df_vcf = pd.concat((df_vcf, mat), axis=1)
    return df_vcf
    

def gt_to_str(h):
    h = h.astype(str)
    h[ np.logical_and(h != '0', h != '1')] = '.'
    return h
    
def _get_gt_and_as(rr, aa, h1, h2):
    h1 = gt_to_str(h1[:, :, np.newaxis])
    h2 = gt_to_str(h2[:, :, np.newaxis])
    hap = np.concatenate((h1, h2), axis=2)
    gt = np.apply_along_axis(_gt, 2, hap)
    
    asc = np.concatenate((rr[:, :, np.newaxis], aa[:, :, np.newaxis])).astype(str)
    as = np.apply_along_axis(_as, 2, asc)
    
def _gt(ii):
    return '|'.join(list(ii))
    
def _as(ii):
    return ','.join(list(ii))

def add_column(df, df_new):
    if df is None:
        return df_new
    else:
        return pd.merge(df, df_new, on='VARIANT_ID', how='outer').fillna(0)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='prep_vcf.py', description='''
        Prepare VCF file for RASQUAL.
    ''')
    parser.add_argument('--asc_prefix', help='''
        Prefix of ASC table
    ''')
    parser.add_argument('--asc_suffix', help='''
        Suffix of ASC table
    ''')
    parser.add_argument('--sample_meta_data', nargs='+', help='''
        Sample meta information. 
        Have the tissue query string in the second place.
    ''')
    parser.add_argument('--genotype_parquet', help='''
        Genotype parquet (with wildcards {chr_num})
    ''')
    parser.add_argument('--chrs', nargs='+', default=None, help='''
        Chromosomes to process.
    ''')
    parser.add_argument('--output', help='''
        output TSV.GZ (with wildcards {chr_num})
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
    
    if args.chrs is None:
        args.chrs = [ i for i in range(1, 23) ]
    
    logging.info('Loading samples')
    df_sample = pd.read_csv(args.sample_meta_data[0], sep='\t')
    df_sample = df_sample[ df_sample['SMTSD'] == args.sample_meta_data[1] ].reset_index(drop=True)
    sample_id_list = df_sample['SAMPID'].unique().to_list()
    indiv_id_list = sample2indiv(sample_id_list)
    if len(sample_id_list) != len(indiv_id_list):
        raise ValueError('Duplicated individual ID.')
    
    logging.info('Going over ASC tables.')
    df_ref = None
    df_alt = None
    for indiv, sample in tqdm(zip(indiv_id_list, sample_id_list)):
        tbl_file = '{}{}{}'.format(args.asc_prefix, indiv, args.asc_suffix)
        df_asc_ = pd.read_csv(tbl_file, sep='\t', compression='gzip')
        df_asc_ = df_asc_[ df_asc_.SAMPLE_ID == sample ].reset_index(drop=True)
        df_r = df_asc_[['VARIANT_ID', 'REF_COUNT']]
        df_r = df_r.rename(columns={indiv: 'REF_COUNT'}, inplace=True)
        df_a = df_asc_[['VARIANT_ID', 'ALT_COUNT']]
        df_a = df_a.rename(columns={indiv: 'ALT_COUNT'}, inplace=True)
        df_ref = add_column(df_ref, df_r)
        df_alt = add_column(df_alt, df_a)
    df_ref = df_ref[ ['VARIANT_ID'] + indiv_id_list ].copy()
    df_alt = df_alt[ ['VARIANT_ID'] + indiv_id_list ].copy()
    
    logging.info('Going over chromosomes')
    for ch in tqdm(args.chrs):
        geno1 = pd.read_parquet(args.genotype_parquet.format(chr_num=ch, hap=1))
        geno2 = pd.read_parquet(args.genotype_parquet.format(chr_num=ch, hap=2))
        geno1 = geno1[ indiv_id_list ].copy()
        geno2 = geno2[ indiv_id_list ].copy()
        df_vcf = construct_vcf(df_ref, df_alt, geno1, geno2)
        df_vcf.to_csv(args.output.format(chr_num=chr), sep='\t', compression='gzip')
    
