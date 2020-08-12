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

snp_list = fread(cmd = paste0('zcat ', opt$snp_list), header = T, sep = '\t')
snp_list$chr = as.character(snp_list$chr)
gene_model = fread(opt$gene_model, header = T, sep = '\t')
gene_model = gene_model[, c('gene_id', 'chromosome', 'strand', 'start', 'end')]
gene_model$chromosome = stringr::str_remove(gene_model$chromosome, 'chr')
tmp = rep(1, nrow(gene_model))
tmp[gene_model$strand == '-'] = -1
gene_model$strand = tmp
colnames(gene_model)[c(2, 4:5)] = c('chr', 'exon_starts', 'exon_ends')
gene_model = gene_model[ gene_model$chr %in% snp_list$chr, ]
gene_model$exon_starts = as.character(gene_model$exon_starts)
gene_model$exon_ends = as.character(gene_model$exon_ends)
gene_model$strand = as.integer(gene_model$strand)
snp_counts = countSnpsOverlapingExons(gene_model, snp_list, cis_window = opt$cis_window_size)
write.table(snp_counts, opt$output, quo = F, row = F, col = T, sep = '\t')

