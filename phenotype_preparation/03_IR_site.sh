#PBS -N Job_Name 
#PBS -q colo-bliu
#PBS -l walltime=36:00:00
#PBS -l select=1:ncpus=1:mem=256GB 
#PBS -o Out_File_Name.out 
#PBS -e Error_File_Name.err 

cd $PBS_O_WORKDIR
python=~/miniconda3/envs/R4.3.1/bin/python
Rscript=~/miniconda3/envs/R4.3.1/bin/Rscript

$Rscript /data/zhangyuntian/project/scSplice/code/new_quantification_method_4_14/03_IR_site.R 
