# microbe pca plot
qiime tools import --type "PCoAResults" \
		   --input-path ../results/microbes-pca.results \
		   --output-path ../results/hfd_output/microbes-pca.qza

qiime emperor plot --i-pcoa ../results/hfd_output/microbes-pca.qza \
		   --m-metadata-file ../../data/HFD/microbe_feature_metadata.txt \
		   --p-ignore-missing-samples \
		   --o-visualization ../results/hfd_output/emperor-microbe-pca.qzv


# metabolite pca plot
# qiime tools import --type "PCoAResults" \
#                    --input-path results/metabolites-pca.results \
#                    --output-path results/metabolites-pca.qza
#
# qiime emperor plot --i-pcoa results/metabolites-pca.qza \
#                    --m-metadata-file data/tsne_metabolite_feature_metadata.txt \
#                    --p-ignore-missing-samples \
#                    --o-visualization results/emperor-metabolite-pca.qzv
#


# biplot
qiime tools import --type "PCoAResults % Properties('biplot')" \
		   --input-path ../results/hfd_output/omics-biplot.results  \
		   --output-path ../results/hfd_output/omics-biplot-results.qza

qiime emperor biplot --i-biplot ../results/hfd_output/omics-biplot-results.qza \
		     --m-sample-metadata-file ../data/HFD/biplot_metabolite_feature_metadata.txt \
		     --m-feature-metadata-file ../data/HFD/microbe_feature_metadata.txt \
		     --p-ignore-missing-samples \
		     --p-number-of-features 100 \
		     --o-visualization ../results/hfd_output/emperor-mouse-omics.qzv
qiime tools view ../results/hfd_output/emperor-mouse-omics.qzv
