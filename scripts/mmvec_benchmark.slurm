#!/bin/bash
#
#SBATCH --job-name=master
#SBATCH --output=stdout.txt
#SBATCH --mail-type=ALL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jmorton@flatironinstitute.org
#
#SBATCH --ntasks=16
#SBATCH --time=30:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH -N 1

echo $PWD
source ~/.bashrc

module load slurm

#module load cuda/8.0.61
#module load cudnn/v7.0-cuda-8.0

source ~/venvs/mmvec-tf/bin/activate
# conda activate mmvec

cd /mnt/home/jmorton/research/multiomics/CF_small_benchmark
directory="../results/benchmark_output/"
# these parameters are for TF
beta1=0.85
beta2=0.9
lr=1e-3
batch=50
dim=3
f="CF_benchmark"


for r in tfrun8 tfrun9 tfrun10
do
	mode='rel'
 	run="rel_${r}_dim${dim}_lr${lr}_beta1${beta1}_beta2${beta2}_batch${batch}"
 	summary="CF_sims/${mode}/mmvec/rel_${r}_results/summary_${run}"
 	model="CF_sims/${mode}/mmvec/rel_model_${run}.txt"
 	ord="CF_sims/${mode}/mmvec/ordination_${run}.txt"
	file1="table.rel.microbes.0.biom"
	file2="table.rel.metabolites.0.biom"
	metadata="metadata.0.txt"
 	sbatch -N 1 -p gpu --gres=gpu:01 --wrap "mmvec paired-omics --microbe-file $relfile1 --metabolite-file $relfile2 --epochs 100 --min-feature-count 10 --learning-rate 1e-3 --beta1 0.85 --beta2 0.9 --batch_size 50 --summary-dir $summary --checkpoint-interval 3600 --summary-interval 600 --arm-the-gpu --training-column Testing --metadata-file $metadata"


	mode='abs'
 	run="rel_${r}_dim${dim}_lr${lr}_beta1${beta1}_beta2${beta2}_batch${batch}"
 	summary="CF_sims/${mode}/mmvec/rel_${r}_results/summary_${run}"
 	model="CF_sims/${mode}/mmvec/rel_model_${run}.txt"
 	ord="CF_sims/${mode}/mmvec/ordination_${run}.txt"
	file1="table.rel.microbes.0.biom"
	file2="table.rel.metabolites.0.biom"
	metadata="metadata.0.txt"

 	run="abs_${r}_dim${dim}_lr${lr}_beta1${beta1}_beta2${beta2}_batch${batch}"
 	summary="abs_${r}_results/summary_${run}"
 	model="abs_model_${run}.txt"
	sbatch -N 1 -p gpu --gres=gpu:01 --wrap "mmvec paired-omics --microbe-file $absfile1 --metabolite-file $absfile2 --epochs 100 --min-feature-count 10 --learning-rate 1e-3 --beta1 0.85 --beta2 0.9 --batch_size 50 --summary-dir $summary --checkpoint-interval 3600 --summary-interval 600 --arm-the-gpu --training-column Testing --metadata-file $metadata"
done
