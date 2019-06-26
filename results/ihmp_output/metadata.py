import pandas as pd
ls taxa_filtered.biom_prior1_songbird_results/
ls taxa_filtered.biom_prior1_songbird_results/summary__dim_prior1_lr_beta10.9_beta20.95_sub_batch10/
taxa_diffs = pd.read_table('taxa_filtered.biom_prior1_songbird_results/summary__dim_prior1_lr_beta10.9_beta20.95_sub_batch10/differentials.tsv', index_col=0)
ls
ls ecs_filtered.biom_prior1_songbird_results/
ms_diffs = pd.read_table('ecs_filtered.biom_prior1_songbird_results/summary__dim_prior1_lr_beta10.9_beta20.95_sub_batch10/differentials.tsv', index_col=0)
ms_diffs
taxa_diffs
taxa_diffs.index.name = 'featureid'
mx_diffs.index.name = 'sampleid'
ms_diffs.index.name = 'sampleid'
ls
ms_diffs.to_csv('metabolite_metadata.txt', sep='\t')
taxa_diffs.to_csv('microbe_metadata.txt', sep='\t')
ls
ms_diffs = pd.read_table('ms_filtered.biom_prior1_songbird_results//summary__dim_prior1_lr_beta10.9_beta20.95_sub_batch10/differentials.tsv', index_col=0)
ms_diffs
ms_diffs.index.name='sampleid'
ms_diffs.to_csv('metabolite_metadata.txt', sep='\t')
