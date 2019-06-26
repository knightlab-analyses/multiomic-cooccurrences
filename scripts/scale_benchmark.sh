# this is for the scale benchmark
exec=../ipynb/src/run_models.py
dir=../results/benchmark_output/scale_benchmark
for mode in abs rel
do
    for tool in pearson spearman
    do
	python $exec run-${tool} \
	       --table1-file $dir/microbe_${mode}.biom \
	       --table2-file $dir/metabolite_${mode}.biom \
	       --output-file $dir/$mode/${mode}_${tool}.txt
    done
done
