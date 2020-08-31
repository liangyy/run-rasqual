library(optparse)

option_list <- list(
    make_option(c("-s", "--snp_list"), type="character", default=NULL,
                help="SNP meta info (TSV.GZ)",
                metavar="character"),
    make_option(c("-g", "--gene_model"), type="character", default=NULL,
                help="gene model",
                metavar="character"),
    make_option(c("-d", "--outdir"), type="character", default=NULL,
                help="directory of output",
                metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="output TXT",
                metavar="character"),
    make_option(c("-c", "--cis_window_size"), type="numeric", default=NULL,
                help="cis window size",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(rasqualTools)
library(data.table)
options(datatable.fread.datatable = F)
library(dplyr)

trim_dot = function(ss) {
  unlist(lapply(strsplit(ss, '\\.'), function(x) { x[1] }))
}

rename_col = function(ss, col_old, col_new) {
  if(col_old %in% colnames(ss)) {
    colnames(ss)[which(col_old == colnames(ss))] = col_new
  } 
  ss
}

load_gene_model = function(filename, snp_chr) {
  gene_model = fread(filename, header = T, sep = '\t')
  if('feature_type' %in% colnames(gene_model)) {
    gene_model = gene_model[ gene_model$feature_type == 'exon', ]
    gene_model$gene_id = trim_dot(gene_model$gene_id)
    gene_model = rename_col(gene_model, 'start_location', 'start')
    gene_model = rename_col(gene_model, 'end_location', 'end')
    gene_model = gene_model %>% group_by(gene_id) %>% summarize(chromosome = chromosome[1], strand = strand[1], 
        start = paste0(start, collapse = ','), end = paste0(end, collapse = ','))
  }
  
  gene_model = gene_model[, c('gene_id', 'chromosome', 'strand', 'start', 'end')]
  gene_model$chromosome = stringr::str_remove(gene_model$chromosome, 'chr')
  tmp = rep(1, nrow(gene_model))
  tmp[gene_model$strand == '-'] = -1
  gene_model$strand = tmp
  colnames(gene_model)[c(2, 4:5)] = c('chr', 'exon_starts', 'exon_ends')
  gene_model = gene_model[ gene_model$chr %in% snp_chr, ]
  gene_model$exon_starts = as.character(gene_model$exon_starts)
  gene_model$exon_ends = as.character(gene_model$exon_ends)
  gene_model$strand = as.integer(gene_model$strand)
  
  gene_model
}

snp_list = fread(cmd = paste0('zcat ', opt$snp_list), header = T, sep = '\t')
snp_list$chr = as.character(snp_list$chr)
gene_model = load_gene_model(opt$gene_model, snp_list$chr)
snp_counts = countSnpsOverlapingExons(gene_model, snp_list, cis_window = opt$cis_window_size)
write.table(snp_counts, opt$output, quo = F, row = F, col = T, sep = '\t')

