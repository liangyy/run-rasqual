library(data.table)
options(datatable.fread.datatable = F)
df = fread("for i in `seq 1 22`; do cat qsub/logs/run-run_kidney_cortex-chr$i.out |grep Resources|tail -n 1 | awk -F'=' '{print $2}'|awk -F',' '{print $1}'; done", header = F)
df2 = fread("for i in `seq 1 22`; do cat qsub/logs/run-run_kidney_cortex-chr$i.out |grep Resources|tail -n 1 | awk -F'=' '{print $4}'|awk -F',' '{print $1}'; done", header = F)
parse_time_in_min = function(s) {
  tmp = strsplit(s, ':')[[1]]
  h = as.numeric(tmp[1]) * 60
  m = as.numeric(tmp[2])
  h + m
}

tot_time = sum(sapply(df$V1, parse_time_in_min))
message('kidney_cortex 100kbp ', 'total cpu time = ', tot_time)
wall_tot_time = sum(sapply(df2$V1, parse_time_in_min))
message('kidney_cortex 100kbp ', 'total wall time = ', wall_tot_time)

df = fread("for i in `seq 1 22`; do cat qsub/logs/run-run_kidney_cortex_1mb-chr$i.out |grep Resources|tail -n 1 | awk -F'=' '{print $2}'|awk -F',' '{print $1}'; done", header = F)
df2 = fread("for i in `seq 1 22`; do cat qsub/logs/run-run_kidney_cortex_1mb-chr$i.out |grep Resources|tail -n 1 | awk -F'=' '{print $4}'|awk -F',' '{print $1}'; done", header = F)
parse_time_in_min = function(s) {
  tmp = strsplit(s, ':')[[1]]
  h = as.numeric(tmp[1]) * 60
  m = as.numeric(tmp[2])
  h + m
}

tot_time = sum(sapply(df$V1, parse_time_in_min))
message('kidney_cortex 1Mbp ', 'total cpu time = ', tot_time)
wall_tot_time = sum(sapply(df2$V1, parse_time_in_min))
message('kidney_cortex 1Mbp ', 'total wall time = ', wall_tot_time)
