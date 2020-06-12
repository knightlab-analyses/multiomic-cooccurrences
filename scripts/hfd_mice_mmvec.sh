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

#source ~/venvs/mmvec/bin/activate
conda activate mmvec

cd /mnt/home/jmorton/research/multiomics/data

beta1=0.9
beta2=0.99
lr=1e-4
batch=10
sub=100
f=HFD


for r in run10 run11 run12
do
    for dim in 3 5
    do
 	run="${r}_${f}_dim${dim}_lr${lr}_beta1${beta1}_beta2${beta2}_sub${sub}_batch${batch}"
 	summary="${r}_${f}_results/summary_${run}"
 	model="${r}_${f}_results/model_${run}.txt"
 	ord="${r}_${f}_results/ordination_${run}.txt"
 	file1="/mnt/home/jmorton/research/multiomics/data/${f}/microbes.biom"
 	file2="/mnt/home/jmorton/research/multiomics/data/${f}/metabolites.biom"
	metadata="/mnt/home/jmorton/research/multiomics/data/${f}/sample-metadata.txt"
 	sbatch -N 1 -p gpu --gres=gpu:01 --wrap "rhapsody mmvec --microbe-file $file1 --metabolite-file $file2 --epochs 10000 --min-feature-count 10 --learning-rate 1e-4 --beta1 0.9 --beta2 0.99 --batch-size 50 --summary-dir $summary --checkpoint-interval 3600 --summary-interval 10 --arm-the-gpu --training-column Testing --metadata-file $metadata"
	done
done
