import pandas as pd
from subprocess import PIPE, Popen

def cmdline(command):
    process = Popen(
        args=command,
        stdout=PIPE,
        shell=True
    )
    return process.communicate()[0].decode('ascii')

def get_tss(strand, s, e):
    tss = 0
    if strand == -1:
        tss = e
    elif strand == 1:
        tss = s
    else:
        raise ValueError('Wrong strand')
    return tss

def gen_gene_index_dict(trc):
    gene_index_dict = {}
    cmd = f'cat {trc} | cut -f 1'
    res = cmdline(cmd)
    res = res.split('\n')
    return { i: idx + 1 for idx, i in enumerate(res) }

def get_sample_size(trc):
    cmd = "cat " + trc + " |head -n 1|awk '{print NF}'"
    res = cmdline(cmd)
    return int(res.split('\n')[0]) - 1  # substract the gene name column

def load_rasqual_output(filename):
    header = [ "Feature_ID" ,"rs_ID" ,"Chromosome" ,"SNP_position" ,"Ref_allele" ,"Alt_allele" ,"Allele_frequency_not_MAF" ,"HWE_Chi-square_statistic" ,"Imputation_quality_score_IA" ,"Log_10_Benjamini-Hochberg_Q-value" ,"Chi_square_statistic_2_x_log_Likelihood_ratio" ,"Effect_size_Pi" ,"Sequencing_mapping_error_rate_Delta" ,"Reference_allele_mapping_bias_Phi" ,"Overdispersion" ,"SNP_ID_within_the_region" ,"No_of_feature_SNPs" ,"No_of_tested_SNPs" ,"No_of_iterations_for_null_hypothesis" ,"No_of_iterations_for_alternative_hypothesis" ,"Random_location_of_ties_tie_lead_SNP;_only_useful_with_-t_option" ,"Log_likelihood_of_the_null_hypothesis" ,"Convergence_status_0success" ,"Squared_correlation_between_prior_and_posterior_genotypes_fSNPs" ,"Squared_correlation_between_prior_and_posterior_genotypes_rSNP" ]
    df = pd.read_csv(filename, header=None, sep='\t')
    df.columns = header
    return df

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='call_rasqual.py', description='''
        Python wrapper to call RASQUAL.
    ''')
    parser.add_argument('--gene_list', nargs='+', help='''
        Gene list along with the column name.
    ''')
    parser.add_argument('--gene_trc', help='''
        It is used for obtaining gene index.
    ''')
    parser.add_argument('--gene_metainfo', help='''
        Meta information of gene needed for RASQUAL call.
        Generated in prep.snmk (using R package rasqualTools).
    ''')
    parser.add_argument('--rasqual_vcf', help='''
        VCF with RASQUAL specification. 
    ''')
    parser.add_argument('--rasqual_exe', help='''
        Executable of RASQUAL.
    ''')
    parser.add_argument('--output', help='''
        Output parquet.
    ''')
    parser.add_argument('--nthread', type=int, default=1, help='''
        Number of threads.
    ''')
    parser.add_argument('--cis_window_size', type=int, default=1000000, help='''
        Cis-window size (surrounding TSS).
    ''')
    parser.add_argument('--trc_bin', help='''
        Total read matrix in BIN.
    ''')
    parser.add_argument('--offset_bin', help='''
        Offset matrix in BIN.
    ''')
    parser.add_argument('--covar_bin', help='''
        Covariate matrix in BIN.
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
    
    df_gene = pd.read_csv(args.gene_metainfo, sep='\t')
    df_gene = df_gene[ df_gene.feature_snp_count > 0 ].reset_index(drop=True)
    gene_list = list(pd.read_csv(args.gene_list[0], sep='\t', compression='gzip')[args.gene_list[1]])
    df_gene = df_gene[ df_gene.gene_id.isin(gene_list) ].reset_index(drop=True)
    
    gene_index_dict = gen_gene_index_dict(args.gene_trc)
    sample_size = get_sample_size(args.gene_trc)
    
    df_all = []
    logging.info('{} genes in total'.format(df_gene.shape[0]))
    for i in range(df_gene.shape[0]):
        # if i > 3:
        #     break
        gene_name = df_gene.iloc[i, 0]
        chr_ = df_gene.iloc[i, 1]
        strand = df_gene.iloc[i, 2]
        gene_start = df_gene.iloc[i, 3]
        gene_end = df_gene.iloc[i, 4]
        nfeature = df_gene.iloc[i, 7]
        ncis = df_gene.iloc[i, 8]
        tss = int(get_tss(strand, start, end))
        cis_start = max(1, tss - args.cis_window_size)
        cis_end = tss + args.cis_window_size
        # take the union of cis-window and gene-body
        region_start = min(gene_start, cis_start)
        region_end = max(gene_end, cis_end)
        call = 'tabix {vcf} {chr_}:{region_start}-{region_end} | {rasqual_exe} \
        -y {trc_bin} \
        -k {offset_bin} \
        -n {sample_size} \
        -j {gene_idx} \
        -l {ncis} \
        -m {nfea} \
        -s {gene_start} \
        -e {gene_end} \
        -f {gene_name} \
        -c {tss} \
        -w {cis_window} \
        --n-threads {nthread} \
        
        -x {covar_bin} > {temp_out}'
        cmd = call.format(
            vcf=args.rasqual_vcf,
            chr_=chr_,
            region_start=region_start,
            region_end=region_end,
            rasqual_exe=args.rasqual_exe,
            trc_bin=args.trc_bin,
            offset_bin=args.offset_bin,
            sample_size=sample_size,
            gene_idx=gene_index_dict[gene_name],
            ncis=ncis,
            nfea=nfeature,
            gene_start=gene_start,
            gene_end=gene_end,
            gene_name=gene_name,
            nthread=args.nthread,
            covar_bin=args.covar_bin,
            tss = tss, 
            cis_window=args.cis_window_size
            temp_out=args.output + f'{gene_name}.temp'
        )
        if not os.path.exists(args.output + f'{gene_name}.temp'):
            logging.info(cmd)
            _ = cmdline(cmd)
        df_this = load_rasqual_output(args.output + f'{gene_name}.temp')
        df_all.append(df_this)
    df_all = pd.concat(df_all, axis=0)
    df_all = df_all.astype({'Chromosome': str}) 
    df_all.to_parquet(args.output)
        
            
        
            
    
