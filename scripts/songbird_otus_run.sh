songbird multinomial \
	 --input-biom ../data/otus_nt.biom \
	 --metadata-file ../data/sample-metadata.txt \
	 --formula "depth" \
	 --epoch 1000 \
	 --batch-size 3 \
	 --summary-dir songbird_cf_depth_otus
# can write "intercept only" as formula="1"
