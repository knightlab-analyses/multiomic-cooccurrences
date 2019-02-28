# microbe pca plot
qiime tools import --type "PCoAResults" \
                   --input-path results/microbes-pca.results \
                   --output-path results/microbes-pca.qza

qiime emperor plot --i-pcoa results/microbes-pca.qza \
                   --m-metadata-file data/tsne_microbe_feature_metadata.txt \
                   --p-ignore-missing-samples \
                   --o-visualization results/emperor-microbe-pca.qzv


# metabolite pca plot
qiime tools import --type "PCoAResults" \
                   --input-path results/metabolites-pca.results \
                   --output-path results/metabolites-pca.qza

qiime emperor plot --i-pcoa results/metabolites-pca.qza \
                   --m-metadata-file data/tsne_metabolite_feature_metadata.txt \
                   --p-ignore-missing-samples \
                   --o-visualization results/emperor-metabolite-pca.qzv



# biplot
# qiime tools import --type "PCoAResults % Properties('biplot')" \
#                    --input-path results/omics-biplot.results  \
#                    --output-path results/omics-biplot-results.qza
# 
# qiime emperor biplot --i-biplot results/omics-biplot-results.qza \
#                      --m-sample-metadata-file data/biplot_metabolite_feature_metadata.txt \
#                      --m-feature-metadata-file data/trimmed_microbe_feature_metadata.txt \
#                      --p-ignore-missing-samples \
#                      --p-number-of-features 100 \
#                      --o-visualization results/emperor-mouse-omics.qzv
# qiime tools view results/emperor-mouse-omics.qzv