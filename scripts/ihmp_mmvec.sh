#!/bin/bash
#
#SBATCH --job-name=ihmp
#SBATCH --output=stdout.txt
#SBATCH --mail-type=ALL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jmorton@flatironinstitute.org
#
#SBATCH --ntasks=16
#SBATCH --time=30:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH -N 1
#SBATCH -p gpu
#SBATCH --gres=gpu:01

echo $PWD
source ~/.bashrc

taxa_filtered.biom
module load slurm
module load cuda/10.0.130_410.48
module load cudnn/v7.6.2-cuda-10.0

#module load cuda/8.0.61
#module load cudnn/v7.0-cuda-8.0

source ~/venvs/mmvec-tf/bin/activate
# conda activate mmvec

cd /mnt/home/jmorton/research/ihmp

# these parameters are for TF
batch=5000
learning_rate=1e-5
outprior=0.1
inprior=0.1
beta1=0.9
beta2=0.95
f=taxa_filtered.biom
# 1. c8_pos
# 2. hilic_pos
# 3. hilic_neg
# 4. c18_neg

for f2 in c8_pos.biom c18_neg.biom
do
     for r in tfrun1 tfrun2 tfrun3
     do
	for dim in 3 5
	do
	    run="${f2}_${r}_dim${dim}_lr${lr}_beta1${beta1}_beta2${beta2}_sub${sub}_batch${batch}"
	    summary="$../results/ihmp_output/mmvec_results/{r}_${f}_results/summary_${run}"
	    model="${r}_results/model_${run}.txt"
	    ord="${r}_results/ordination_${run}.txt"
	    file1="../results/ihmp_output/$f"
	    file2="../results/ihmp_output/$f2"
	    sbatch -N 1 -p gpu --gres=gpu:02 --wrap "mmvec paired-omics --microbe-file $file1 --metabolite-file $file2 --epochs 3000 --latent-dim ${dim} --min-feature-count 10 --learning-rate 1e-5 --input_prior 1 --output_prior 1 --beta1 0.9 --beta2 0.95 --batch-size $batch  --summary-dir $summary --checkpoint-interval 3600 --summary-interval 1200 --arm-the-gpu --training-column Testing --metadata-file sample-metadata.txt"

	done
    done
done
nppp
