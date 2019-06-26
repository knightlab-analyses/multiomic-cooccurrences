batch=5
epochs=5000
d=0.5
for f in filtered_otus_nt.biom filtered_lcms_nt.biom
do
	run="${f}_dim${dim}_prior${d}_lr${lr}_beta1${beta1}_beta2${beta2}_sub${sub}_batch${batch}_depth"
	summary="../results/cf_output/${f}_prior${d}_songbird_results/summary_${run}"
	model="../../results/cf_output/${f}_prior${d}_songbird_results/model_${run}.txt"
	file="../../data/CF/$f"
	metadata="../../data/CF/sample-metadata.txt"
	echo $file
 	songbird multinomial --input-biom ${file} --metadata-file $metadata --formula "depth"  --epochs $epochs  --learning-rate 1e-3 --differential-prior $d --batch-size $batch  --summary-dir $summary --checkpoint-interval 3600 --summary-interval 1 --training-column Testing --min-sample-count 1000 --min-feature-count 10

done
