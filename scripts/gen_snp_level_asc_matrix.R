library(optparse)

option_list <- list(
    make_option(c("-a", "--asc"), type="character", default=NULL,
                help="ASC SNP level table (by individual)",
                metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="output",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)



library(data.table)
options(datatable.fread.datatable = F)
library(dplyr)


wasp = fread(cmd = paste0('zcat ', opt$asc), header = T, sep = '\t')
wasp = wasp %>% filter(!is.na(GENE_ID), GENOTYPE %in% c('GT;0|1', 'GT;1|0'))
wasp = wasp[, c('VARIANT_ID', 'SAMPLE_ID', 'REF_COUNT', 'ALT_COUNT')]
gz1 <- gzfile(opt$output, "w")
write.table(wasp, gz1, col = T, row = F, quo = F, sep = '\t')
close(gz1)
