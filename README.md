
# RA-Endotypes

This repository is for the project **"Multi-modal integration of protein interactomes with genomic and molecular data discovers distinct RA endotypes"**.

## Overview

Rheumatoid arthritis (RA) is a complex autoimmune disease characterized by clinical and molecular heterogeneity, especially in the presence of anti-cyclic citrullinated peptide antibodies (CCP). CCP positivity in RA patients is associated with more severe disease progression and distinct treatment responses compared to CCP-negative patients. Although cellular and molecular differences between RA subtypes have been studied, genetic differences at a network scale remain largely unexplored.

Here, we use a novel multi-scale framework that integrates a **network-based genome-wide association study (GWAS)** with **functional genomic data** to identify network modules distinguishing CCP+ and CCP– RA. Using the RACER (Rheumatoid Arthritis Comparative Effectiveness Research) cohort, comprising 555 CCP+/RF+ and 384 CCP-/RF+ RA patients, we identified 14 gene modules that explain genetic differences between CCP+/RF+ and CCP-/RF+ RA. These findings were validated through heritability partitioning and multivariate expression analyses, uncovering novel genetic loci associated with RA heterogeneity and supporting the development of more personalized therapeutic strategies.

![GitHub Front Image](https://github.com/user-attachments/assets/d1a68c65-4379-4726-b8c3-884a0cfaba9c)

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Reproducibility](#reproducibility)
- [Contributing](#contributing)
- [License](#license)

## Installation

To set up the environment and run the code, follow these steps:

1. **Clone the repository:**
   ```bash
   git clone https://github.com/jishnu-lab/RAEndo.git
   cd RAEndo
   ```
2. **Install the required dependencies:**
   
| Methods                          | Instruction                                                 |
|----------------------------------------|-------------------------------------------------------------|
| **PrediXcan** | https://github.com/hakyimlab/PrediXcan           |
| **heritability**| https://dougspeed.com/downloads2/              |
| **hotnet** | https://raphael-group.github.io/projects/hotnet2/       |

## Usage

### Running the Analysis

This analysis involves running multiple scripts in sequence. There is no main script, so each script should be executed individually. Follow the order outlined below for the complete analysis:

1. **Heritability**: 
   - This script estimates heritability across RA subtypes using genomic data.
   - Navigate to the `heritability` directory:
     ```bash
     cd Code/heritability
     ```
   - Run the script:
     ```bash
     bash make_grms.sh
     bash make_estimation.sh
     
     ```

2. **HotNet2**: 
   - This script identifies subnetworks associated with RA subtypes using the HotNet2 algorithm.
   - Navigate to the `hotnet` directory:
     ```bash
     cd Code/hotnet
     ```
   - Run the script:
     ```bash
     bash Run_hotnet.sh
     Rscript process_hotnet_updated.R
     ```

3. **Bulk RNA-Seq Analysis**: 
   - This script analyzes bulk RNA-seq data to identify gene expression patterns in RA subtypes.
   - Navigate to the `bulk-rna` directory:
     ```bash
     cd Code/bulk-rna
     ```
   - Run the script:
     ```bash
     Rscript HLA_beta_feat_perm.R
     ```

4. **PrediXcan**: 
   - This script performs gene expression prediction from genotypes using PrediXcan.
   - Navigate to the `PrediXcan` directory:
     ```bash
     cd Code/PrediXcan
     ```
   - Run the script:
     ```bash
     Rscript imputation.R
     ```

Ensure you execute each script in the order listed above to maintain consistency and reproducibility of results.

- For individual analyses, see the [Reproducibility](#reproducibility) table below.

## Reproducibility

The table below provides links to the Jupyter notebooks used in the analysis. Each notebook corresponds to a specific step in the project, allowing you to reproduce the results easily.


| Notebook Name                          | Description                                                 |
|----------------------------------------|-------------------------------------------------------------|
| [**PrediXcan**](https://github.com/jishnu-lab/RAEndo/tree/main/Code/PrediXcan) | Performs PrediXcan analysis to predict gene expression based on genotypes.            |
| [**bulk-rna**](https://github.com/jishnu-lab/RAEndo/tree/main/Code/bulk-rna) | Analyzes bulk RNA-seq data to identify expression patterns in RA subtypes.            |
| [**heritability**](https://github.com/jishnu-lab/RAEndo/tree/main/Code/heritability) | Computes heritability estimates using genomic data across RA subtypes.               |
| [**hotnet**](https://github.com/jishnu-lab/RAEndo/tree/main/Code/hotnet) | Applies HotNet2 algorithm to detect subnetworks associated with RA subtypes.         |

## Contributing

We welcome contributions to improve this project! Please open an issue or submit a pull request with your proposed changes.

## License

This project is licensed under the [MIT License](LICENSE).
