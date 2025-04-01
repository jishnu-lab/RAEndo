plink --bfile dataset_43345835_merged-maf-0.01-geno-0.1 --hwe 0.001 --make-bed --out dataset_43345835_merged-maf-0.01-geno-0.1_hwe_0.001

plink --bfile dataset_43345835_merged_maf_05_1_hwe_0.001 --pheno updated_file.pheno --covar gender.cov --logistic --hide-covar --out dataset_43345835_merged_maf_05_1_hwe_0.001 --allow-no-sex

plink --maf 0.05 --geno 0.1 --make-bed --bfile dataset_43345835_merged --out ./plink_files/dataset_43345835_merged_maf_05_
