library(data.table)
options(datatable.fread.datatable = F)
df = fread("for i in `seq 1 22`; do cat qsub/logs/run-run_kidney_cortex-chr$i.out |grep Resources|tail -n 1 | awk -F'=' '{print $2}'|awk -F',' '{print $1}'; done", header = F)

parse_time_in_min = function(s) {
  tmp = strsplit(s, ':')[[1]]
  h = as.numeric(tmp[1]) * 60
  m = as.numeric(tmp[2])
  h + m
}

tot_time = sum(sapply(df$V1, parse_time_in_min))
message(tot_time)

