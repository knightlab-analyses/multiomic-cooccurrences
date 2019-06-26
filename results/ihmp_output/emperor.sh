
for mode in c18_neg c8_pos hilic_pos hilic_neg
do
    qiime tools import \
	  --input-path ${mode}.biom.ordination.txt \
	  --output-path ${mode}_ordination.qza \
	  --type "PCoAResults % Properties('biplot')"

    qiime emperor biplot \
	  --i-biplot ${mode}_ordination.qza \
	  --m-sample-metadata-file ${mode}.biom.metadata.txt \
	  --m-feature-metadata-file microbe_metadata.txt \
	  --p-ignore-missing-samples \
	  --p-number-of-features 50 \
	  --o-visualization ${mode}_emperor.qzv

done

