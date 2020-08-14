library(optparse)

option_list <- list(
    make_option(c("-t", "--trc"), type="character", default=NULL,
                help="TRC bed file",
                metavar="character"),
    make_option(c("-c", "--covar"), type="character", default=NULL,
                help="covariate",
                metavar="character"),
    make_option(c("-i", "--indiv_list"), type="character", default=NULL,
                help="the list of individual id. it defines the order. (read header of vcf)",
                metavar="character"),
    make_option(c("-d", "--outdir"), type="character", default=NULL,
                help="directory of output",
                metavar="character"),
    make_option(c("-o", "--output_trc_prefix"), type="character", default=NULL,
                help="prefix of TRC matrix BIN",
                metavar="character"),
    make_option(c("-w", "--output_sf_prefix"), type="character", default=NULL,
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

is_match = function(str, pat) {
  tmp = stringr::str_match(str, pat)
  return(!is.na(tmp))
}

load_header = function(vcf) {
  e = fread(cmd = paste0('zcat < ', vcf, ' | head -n 2')) 
  colnames(e)[-(1:9)]
}



indiv_list = load_header(opt$indiv_list)

trc = fread(cmd = paste0('zcat ', opt$trc), header = T, sep = '\t')
trc_mat = trc[, c(-1, -2, -3, -4)]
trc_mat = trc_mat[, indiv_list]
trc_mat = as.matrix(trc_mat)
rownames(trc_mat) = trc[, 4]
# remove constant gene
tmp = apply(trc_mat, 1, sd)
trc_mat = trc_mat[ tmp != 0, ]
name_tag = strsplit(opt$output_trc_prefix, '\\.')[[1]]
olist = list()
olist[[name_tag[1]]] = trc_mat
saveRasqualMatrices(olist, opt$outdir, file_suffix = name_tag[2])

size_factor = rasqualCalculateSampleOffsets(trc_mat, gc_correct = FALSE)
name_tag = strsplit(opt$output_sf_prefix, '\\.')[[1]]
olist = list()
olist[[name_tag[1]]] = size_factor
saveRasqualMatrices(olist, opt$outdir, file_suffix = name_tag[2])

covar_pc = rasqualTools:::rasqualMakeCovariates(trc_mat, size_factor)
covar_in = fread(cmd = paste0('zcat ', opt$covar))
is_peer = is_match(covar_in[, 1], 'InferredCov')
covar_in = covar_in[ !is_peer, ]
covar_mat = covar_in[, -1]
covar_mat = t(as.matrix(covar_mat))
covar_mat = covar_mat[ colnames(trc_mat), ]
covar_mat = cbind(covar_mat, covar_pc)
name_tag = strsplit(opt$output_covar_prefix, '\\.')[[1]]
olist = list()
olist[[name_tag[1]]] = covar_mat
saveRasqualMatrices(olist, opt$outdir, file_suffix = name_tag[2])

