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

## Detailed Description of Functions

### `arg_parser()`

Parses command-line arguments required for running the script.

### `print_header()`

Prints the script header with author and affiliation information.

### `reverse_complement(seq)`

Returns the reverse complement of a given DNA sequence.

### `check_dependencies()`

Checks for the presence of `fastp` and `weblogo` dependencies.

### `print_status(msg)`

Prints status messages with timestamps.

### `open_file(file)`

Opens a file, handling both gzipped and regular text files.

### `run_PAM_count(FastqFile, upOutputFile, downOutputFile, TargetSeq, UpstreamLen, DownstreamLen, SeqDirection)`

Counts the PAM sequences in a Fastq file around a target sequence.

### `Generate_heatmap(forMerCount, revMerCount, totalMerCount, heatpng, flank, flanklen)`

Generates a heatmap for PAM sequences.

### `Draw_PAM_Heatmap(pamcount, heatpng, flank, flanklen)`

Draws a heatmap for PAM counts.

### `get_file_name_dict(ExpFileList, ConFileList, PAMCountFileTxt)`

Creates a dictionary of file names for experiment and control samples.

### `Merge_PAM_count(PAMCountList, MergedFile, file_name_dict, upstreamLen)`

Merges PAM count files into a single file.

### `read_and_prepare_data(rawfile)`

Reads and preprocesses PAM count data, generating statistical features.

### `train_model(X_train, y_train)`

Trains a Random Forest classifier and performs grid search for hyperparameter optimization.

### `evaluate_model(clf, X_test, y_test, eval_roc_png, eval_prc_png)`

Evaluates the trained model and generates ROC and PRC curves.

### `save_classification_report(clf, X_test, y_test, report_file)`

Saves the classification report of the trained model.

### `generate_logo_data(data, PAMlen)`

Generates data for creating PAM logos.

### `save_logo_data(logo_data, resfile)`

Saves the data required for creating PAM logos.

### `Combine_PAM_score(upstreamPAMscore, downstreamPAMscore, compinedPAMscore, UpstreamLen, DownstreamLen)`

Combines upstream and downstream PAM scores for unified visualization.

### `create_pam_logo_combined(compinedPAMscore, combinedPAMlogo, UpstreamLen, DownstreamLen)`

Creates a combined PAM logo using `weblogo`.

### `create_pam_logo(PAMscoreFile, PAMlogoFile, directions, PAMlen, overall_confidence_score)`

Creates a PAM logo for upstream or downstream sequences.

### `visualize_fastp_results(fastp_dir, output_dir)`

Visualizes the results of `fastp` filtering.

### `Process_file(PAMCountFile, PAMscoreFile, PAMlogoFile, evalROC, evalPRC, clfReportFile, PAMlen, directions, tempfile)`

Processes PAM count files, generates logos, and evaluates models.


## Acknowledgements

This tool was developed by Chen QI, Baitao LI, and collaborators from BNU-HKBU United International College, University of Chinese Academy of Sciences, and BGI Research.

For any questions or support, please contact [qichen@uic.edu.cn].

```
