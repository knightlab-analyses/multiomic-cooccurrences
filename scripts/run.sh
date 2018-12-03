
# parameters that you will need to tune sorted in terms of priority
CURPATH=../
microbe_file=../data/otus_nt.biom
metabolite_file=../data/lcms_nt.biom
latent_dim=3
learning_rate=1e-5
epochs=1000
outprior=1
inprior=1
beta1=0.85
beta2=0.90
batch_size=1000

RESULTS_DIR=$CURPATH/results/latent_dim_${latent_dim}_in_${inprior}_out_${outprior}_beta1_${beta1}_beta2_${beta2}
minstrel autoencoder \
      --microbe-file $microbe_file \
      --metabolite-file $metabolite_file \
      --min-feature-count 10 \
      --num-testing-examples 10 \
      --learning-rate $learning_rate \
      --input_prior $inprior \
      --output_prior $outprior \
      --latent_dim $latent_dim \
      --beta1 $beta1 \
      --beta2 $beta2 \
      --summary-dir $RESULTS_DIR \
      --batch_size $batch_size \
      --epochs $epochs \
      --threads 1 \
      --top-k 10 \
      --summary-interval 10 \
      --checkpoint-interval 3600
