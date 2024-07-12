# CUTPLASMID: Cleavage Under Targeting and PLAsmid Sequencing for Motif IDentification

## Introduction

**CUTPLASMID** is a comprehensive tool designed for the identification and analysis of PAM (Protospacer Adjacent Motif) sequences in genome data, crucial for CRISPR-Cas systems. This script processes genomic data, identifies PAM sequences, and provides detailed statistical and visual analyses, including heatmaps and machine learning-based classification reports.

## Features

- **Data Preprocessing**: Supports raw data filtering using `fastp`.
- **PAM Sequence Counting**: Counts upstream and downstream PAM sequences around a target sequence.
- **Heatmap Generation**: Creates heatmaps to visualize PAM sequence distributions.
- **Machine Learning**: Trains and evaluates models to classify significant PAM sequences.
- **Visualization**: Generates various plots and logos to represent PAM sequence statistics and distributions.

## Dependencies

Ensure the following dependencies are installed:

- Python 3.6+
- pandas
- numpy
- argparse
- datetime
- itertools
- gzip
- json
- seaborn
- Bio (Biopython)
- scipy
- sklearn
- matplotlib
- fastp
- weblogo

Install dependencies using pip:
```sh
pip install pandas numpy argparse datetime biopython scipy scikit-learn matplotlib seaborn
```

## Installation

1. Clone the repository:
```sh
git clone https://github.com/yourusername/CUTPLASMID.git
cd CUTPLASMID
```

2. Ensure you have `fastp` and `weblogo` installed and available in your PATH.

## Usage

### Command-Line Arguments

| Argument         | Description                                                                 |
|------------------|-----------------------------------------------------------------------------|
| `-e, --experiment` | The experiment data file, separated by commas (required)                     |
| `-c, --control`   | The control data file, separated by commas (required)                        |
| `-o, --output`    | The output directory (required)                                              |
| `-t, --target`    | The target sequence (required)                                               |
| `-u, --upstream`  | The upstream PAM library length (default: 4)                                  |
| `-d, --downstream`| The downstream PAM library length (default: 4)                                |
| `-f, --fastp`     | Use fastp to filter the raw data (optional)                                   |
| `-sd, --seqDirection` | The sequence direction of the PAM (`forward` or `reverse`, default: `forward`) |
| `-hm, --heatmap`  | Draw the heatmap of the PAM count (optional)                                  |
| `-kt, --keeptmp`  | Keep the temporary files (optional)                                           |

### Example

```sh
python cutplasmid.py -e experiment1.fastq,experiment2.fastq -c control1.fastq,control2.fastq -o ./output -t AGGCTAGC -u 4 -d 4 -f -hm
```

## Acknowledgements

This tool was developed by Chen QI, Baitao LI, and collaborators from BNU-HKBU United International College, University of Chinese Academy of Sciences, and BGI Research.

For any questions or support, please contact [qichen@uic.edu.cn].
