

#PBS -N Job_Name 
#PBS -q colo-bliu
#PBS -l walltime=36:00:00
#PBS -l select=1:ncpus=1:mem=256GB 
#PBS -o Out_File_Name.out 
#PBS -e Error_File_Name.err 

cd $PBS_O_WORKDIR
python=~/miniconda3/envs/R3.6.0/bin/python
Rscript=~/miniconda3/envs/R3.6.0/bin/Rscript

$python 03_site_identify.py \
  --cluster_file /data/zhangyuntian/simulation/PSI_centric/test/leafcutter_refined \
  --run_dir /data/zhangyuntian/project/scSplice/code/new_quantification_method_4_14/result \
  --prefix test

$python 04_compute_PSI.py \
  --sample /data/zhangyuntian/simulation/PSI_centric/test/leafcutter_perind_numers.constcounts.gz \
  --splice_site /data/zhangyuntian/project/scSplice/code/new_quantification_method_4_14/result/test \
  --run_dir /data/zhangyuntian/project/scSplice/code/new_quantification_method_4_14/result \
  --prefix test
  

$Rscript 05_filter_phenotype.R \
 --input_dir /data/zhangyuntian/project/scSplice/code/new_quantification_method_4_14/result \
 --input_prefix  test \
 --output_dir /data/zhangyuntian/project/scSplice/code/new_quantification_method_4_14/result \
 --output_prefix test_reform.txt \
 -c /data/zhangyuntian/simulation/PSI_centric/test/leafcutter_perind.constcounts.gz \
 -p 0.6 \
 -s 0.1

for i in `seq 2 4`;do
 Rscript 06_DS.R \
  -i /data/zhangyuntian/simulation/intron_centric/group_${i}.txt \
  --input_dir /data/zhangyuntian/project/scSplice/code/new_quantification_method_4_14/result \
  --input_prefix test_reform.txt \
  --output_dir /data/zhangyuntian/project/scSplice/code/new_quantification_method_4_14/result \
  --output_prefix test_G${i}_sGene
done
