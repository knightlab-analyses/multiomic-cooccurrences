batch=5
epochs=5000
d=0.5
for f in filtered_otus_nt.biom
do
	run="${f}_dim${dim}_prior${d}_lr${lr}_beta1${beta1}_beta2${beta2}_sub${sub}_batch${batch}_id_depth"
	summary="../results/cf_output/${f}_prior${d}_songbird_results/summary_${run}"
	model="../results/cf_output/${f}_prior${d}_songbird_results/model_${run}.txt"
	file="../data/CF/$f"
	metadata="../data/CF/sample-metadata.txt"
	#formula="diagnosis+sex+interval_days"
	echo $file
 	songbird multinomial --input-biom ${file} --metadata-file $metadata --formula "C(host_subject_id) + depth"  --epochs $epochs  --learning-rate 1e-3 --differential-prior $d --batch-size $batch  --summary-dir $summary --checkpoint-interval 3600 --summary-interval 1 --training-column Testing --min-sample-count 1000 --min-feature-count 10

done

deicode --in-biom ../data/CF/filtered_otus_nt.biom --output-dir ../results/cf_output/deicode_microbes --rank 3 --iterations 10
deicode --in-biom ../data/CF/filtered_lcms_nt.biom --output-dir ../results/cf_output/deicode_metabolites --rank 3 --iterations 10
