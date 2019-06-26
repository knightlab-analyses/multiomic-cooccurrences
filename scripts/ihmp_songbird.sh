#!/bin/bash
#
#SBATCH --job-name=micetrain
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


# conda activate mmvec

#cd ../data/ihmp/

# these parameters are for TF
batch=10
learning_rate=1e-3
beta1=0.9
beta2=0.95

# f=ecs_filtered.biom
# d=1
# run="dim${dim}_prior${d}_lr${lr}_beta1${beta1}_beta2${beta2}_sub${sub}_batch${batch}"
# summary="${f}_prior${d}_songbird_results/summary_${run}"
# model="${f}_prior${d}_songbird_results/model_${run}.txt"
# file="$f"

#songbird multinomial --input-biom ${file} --metadata-file sample-metadata.txt --formula "diagnosis" --learning-rate 1e-3 --differential-prior $d --batch-size $batch --min-sample-count 1 --min-feature-count 1 --summary-dir $summary --checkpoint-interval 3600 --summary-interval 1200 --training-column "Testing"

#for f in ecs_filtered.biom taxa_filtered.biom ms_filtered.biom
for f in taxa_filtered.biom hilic_pos.biom hilic_neg.biom c8_pos.biom c18_neg.biom
do
    for d in 0.1 1
    do
	run="${r}_dim${dim}_prior${d}_lr${lr}_beta1${beta1}_beta2${beta2}_sub${sub}_batch${batch}"
	summary="../results/ihmp_output/${f}_prior${d}_songbird_results/summary_${run}"
	model="../results/ihmp_output/${f}_prior${d}_songbird_results/model_${run}.txt"
	file="../data/ihmp/$f"
	metadata="../data/ihmp/sample-metadata.txt"
	#formula="diagnosis+sex+interval_days"
	echo $file

 	songbird multinomial --input-biom ${file} --metadata-file $metadata --formula "diagnosis" --learning-rate 1e-3 --differential-prior $d --batch-size $batch  --summary-dir $summary --checkpoint-interval 3600 --summary-interval 1200 --training-column Testing --min-sample-count 1 --min-feature-count 1

    done
done
