data_dir=data
picrust2_pipeline.py \
    -s $data_dir/reference-hit.seqs.fa \
    -i $data_dir/reference-hit.biom \
    -o results/picrust_out_pipeline \
    --threads 8
