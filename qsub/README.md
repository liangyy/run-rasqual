* Step1: `qsub -v MYCONF=vcf_whole_blood run_prep_vcf.qsub` 
* Step2: `qsub -v MYCONF=prep_whole_blood,CHR_NUM=22 run_prep.qsub` (process chr22)
* Step3: `qsub -v MYCONF=run_whole_blood,CHR_NUM=22 run.qsub` (process chr22)
