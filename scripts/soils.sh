# this is for the scale benchmark
exec=../ipynb/src/run_models.py
dir=../data/soils
output=../results/soil_output

for tool in pearson spearman
do
    python $exec run-${tool} \
	   --table1-file $dir/microbes.biom \
	   --table2-file $dir/metabolites.biom \
	   --output-file $output/soils_${tool}.txt
done
