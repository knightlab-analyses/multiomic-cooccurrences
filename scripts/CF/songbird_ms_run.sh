songbird multinomial \
	 --input-biom ../data/lcms_nt.biom \
	 --metadata-file ../data/sample-metadata.txt \
	 --formula "depth" \
	 --epoch 1000 \
	 --batch-size 3 \	 
	 --summary-dir songbird_cf_depth_lcms
