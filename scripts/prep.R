library(optparse)

option_list <- list(
    make_option(c("-t", "--trc"), type="character", default=NULL,
                help="TRC bed file",
                metavar="character"),
    make_option(c("-c", "--covar"), type="character", default=NULL,
                help="covariate",
                metavar="character"),
    make_option(c("-d", "--outdir"), type="character", default=NULL,
                help="directory of output",
                metavar="character"),
    make_option(c("-o", "--output_trc_prefix"), type="character", default=NULL,
                help="prefix of TRC matrix BIN",
                metavar="character"),
    make_option(c("-u", "--output_sf_prefix"), type="character", default=NULL,
                help="prefix of score factor BIN",
                metavar="character"),
    make_option(c("-u", "--output_covar_prefix"), type="character", default=NULL,
                help="prefix of covariate BIN",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(rasqualTools)
library(data.table)
options(datatable.fread.datatable = F)


trc = fread(cmd = paste0('zcat ', opt$trc), header = T, sep = '\t')
trc_mat = trc[, c(-1, -2, -3)]
name_tag = strsplit(opt$output_trc_prefix, '\\.')[[1]]
olist = list()
oilst[[name_tag[1]]] = trc_mat
saveRasqualMatrices(oilst, opt$outdir, file_suffix = name_tag[2])

size_factor = rasqualCalculateSampleOffsets(trc_mat, gc_correct = FALSE)
name_tag = strsplit(opt$output_sf_prefix, '\\.')[[1]]
olist = list()
oilst[[name_tag[1]]] = size_factor
saveRasqualMatrices(olist, opt$outdir, file_suffix = name_tag[2])

covar_pc = rasqualTools:::rasqualMakeCovariates(trc_mat, size_factor)
covar_in = fread(cmd = paste0('zcat ', opt$))
covar_mat = covar[, -1]
covar_mat = covar_mat[, colnames(trc_mat)]
covar_mat = rbind(covar_mat, covar_pc)
name_tag = strsplit(opt$output_covar_prefix, '\\.')[[1]]
olist = list()
oilst[[name_tag[1]]] = covar_mat
saveRasqualMatrices(oilst, opt$outdir, file_suffix = name_tag[2])
