ldak5.2.linux --cut-genes cut_genes  --bfile racer-mega-imputed-f-matched-v3-QC --genefile LDAK_annot.txt
ldak5.2.linux  --calc-genes-reml cut_genes --pheno ../racer-mega-imputed-f-matched-v3-QC.fam  --mpheno 4  --bfile  ../racer-mega-imputed-f-matched-v3-QC  --ignore-weights YES --power -.25
--covar racer-mega-imputed-f-matched-v3-QC-5PC-gender.cov --max-threads 30
ldak5.2.linux --join-genes-reml cut_genes
