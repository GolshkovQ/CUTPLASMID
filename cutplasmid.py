import sys, os
import pandas as pd
import numpy as np
import argparse
import datetime
import itertools
import gzip
import json
import seaborn as sns
from Bio import SeqIO
from scipy.stats import zscore, kruskal, ttest_ind
from sklearn.preprocessing import RobustScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_predict, GridSearchCV, train_test_split
from sklearn.metrics import roc_curve, auc, precision_recall_curve, confusion_matrix, classification_report
import matplotlib.pyplot as plt

def arg_parser():

    ap = argparse.ArgumentParser(description='CUTPLASMID: Cleavage Under Targeting and PLAsmid Sequencing for Motif IDentification')

    ap.add_argument("-e",
                    "--experiment", 
                    required=True, 
                    help="The experiment data file, separated by comma")
    

    ap.add_argument("-c",
                    "--control", 
                    required=True, 
                    help="The control data file, separated by comma")
    
    ap.add_argument("-o",
                    "--output", 
                    required=True, 
                    help="The output directory")
    
    ap.add_argument("-t",
                    "--target", 
                    required=True, 
                    help="The target sequence")
    
    ap.add_argument("-u",
                    "--upstream", 
                    help="The upstream PAM library length", 
                    default=4)
    
    ap.add_argument("-d",
                    "--downstream", 
                    help="The downstream PAM library length", 
                    default=4)
    
    ap.add_argument("-f",
                    "--fastp",
                    help="Use fastp to filter the raw data, if yes, use default parameter. If want modify fastp please edit this script", 
                    action="store_true")
    
    ap.add_argument("-sd",
                    "--seqDirection",
                    choices=['forward', 'reverse'],
                    help="The sequence direction of the PAM, forward means PAM did not need to RC, reverse means PAM need to RC",
                    default='forward')
    
    ap.add_argument("-hm",
                    "--heatmap",
                    action="store_true",
                    help="Draw the heatmap of the PAM count")
    
    ap.add_argument("-kt",
                    "--keeptmp",
                    action="store_true",
                    help="Keep the temporary files")

    opts = ap.parse_args()

    return opts

def print_header():
    print(
        '''
        ================================================================                                                           
                                _                            
              __       -/-,_    // __,   ,    ,____,   .  __/ 
            _(_,__(_/__/__/_)__(/_(_/(__/_)__/ / / (__/__(_/(_
                         /                                    
                        /                                                                                                                                                                               
        CUTPLASMID
        Cleavage Under Targeting and PLAsmid Sequencing for Motif IDentification
        Author: Chen QI, Baitao LI
        BNU-HKBU United International College, Zhuhai, China
        University of Chinese Academy of Sciences, Beijing, China
        BGI Research, Shenzhen, China
        ================================================================
    '''
    )
    return True

def reverse_complement(seq):
    trantab = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(trantab)[::-1]

def check_dependencies():
    ### Fastp dependency
    if os.system("which fastp") != 0:
        print("fastp is not found, please install fastp first")
        sys.exit(1)
    
    ### WebLOGO dependency
    if os.system("which weblogo") != 0:
        print("weblogo is not found, please install weblogo first")
        sys.exit(1)

    return True

def print_status(msg):
    bluecolor = "\033[94m"
    resetcolor = "\033[0m" 

    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    message = bluecolor + "[" + current_time + "] " + resetcolor + msg
    print(message)
    return True

def open_file(file):
    if file.endswith('.gz'):
        return gzip.open(file, "rt")
    else:
        return open(file, "rt")

def run_PAM_count(FastqFile, upOutputFile, downOutputFile, TargetSeq, UpstreamLen, DownstreamLen, SeqDirection):
    upPAMseqs = [''.join(i) for i in itertools.product("ACGT", repeat=UpstreamLen)]
    downPAMseqs = [''.join(i) for i in itertools.product("ACGT", repeat=DownstreamLen)]

    fastqHandle = open_file(FastqFile)
    UpstreamCounterDict = {}
    DownstreamCounterDict = {}
    for PAMseq in upPAMseqs:
        UpstreamCounterDict[PAMseq] = 0
    for PAMseq in downPAMseqs:
        DownstreamCounterDict[PAMseq] = 0
    
    totalCount = 0

    for records in SeqIO.parse(fastqHandle, "fastq"):
        if str(TargetSeq) in str(records.seq):
            targetIndex = str(records.seq).index(str(TargetSeq))
            tempupflankSeq = str(records.seq)[targetIndex-UpstreamLen:targetIndex]
            tempdownflankSeq = str(records.seq)[targetIndex+len(str(TargetSeq)):targetIndex+len(str(TargetSeq))+DownstreamLen]
            if SeqDirection == "forward":
                upflankSeq = tempupflankSeq
                downflankSeq = tempdownflankSeq
            else:
                upflankSeq = reverse_complement(tempdownflankSeq)
                downflankSeq = reverse_complement(tempupflankSeq)
            if len(upflankSeq) == UpstreamLen and len(downflankSeq) == DownstreamLen and not "N" in upflankSeq and not "N" in downflankSeq:
                UpstreamCounterDict[upflankSeq] += 1
                DownstreamCounterDict[downflankSeq] += 1
        totalCount += 1
        if totalCount % 1000000 == 0:
            print_status("Processed "+str(totalCount)+" reads, fastq file: "+FastqFile)

    with open(upOutputFile, "w") as upoutHandle, open(downOutputFile, "w") as downoutHandle:
        upoutHandle.write("PAMseq\tUpstreamCount\n")
        for upPAMseq in upPAMseqs:
            upoutHandle.write(upPAMseq+"\t"+str(UpstreamCounterDict[upPAMseq])+"\n")
        downoutHandle.write("PAMseq\tDownstreamCount\n")
        for downPAMseq in downPAMseqs:
            downoutHandle.write(downPAMseq+"\t"+str(DownstreamCounterDict[downPAMseq])+"\n")
    
    return True

def Generate_heatmap(forMerCount, revMerCount, totalMerCount, heatpng, flank, flanklen):
    bases = ['A', 'C', 'G', 'T']
    forKmer = [''.join(i) for i in itertools.product(bases, repeat=len(list(forMerCount.keys())[0]))]
    revKmer = [''.join(i) for i in itertools.product(bases, repeat=len(list(revMerCount.keys())[0]))]

    heatmap_data = pd.DataFrame(0, index=forKmer, columns=revKmer).fillna(0)
    annot_data = pd.DataFrame(0, index=forKmer, columns=revKmer, dtype=str).fillna("")

    for forMer in forMerCount.keys():
        for revMer in revMerCount.keys():
            totalMer = forMer+revMer
            heatmap_data.at[forMer, revMer] = totalMerCount.get(totalMer, 0)
            freq = totalMerCount.get(totalMer, 0)
            annot_data.at[forMer, revMer] = f"{totalMer}\n{freq}"
    
    plt.figure(figsize=(15,12))
    sns.heatmap(heatmap_data, annot=annot_data, fmt='', cmap='Greys', cbar=False)
    plt.xlabel('Prev-mer')
    plt.ylabel('Suff-mer')
    plt.title(flank)
    plt.savefig(heatpng)
    plt.close()

    return True

def Draw_PAM_Heatmap(pamcount, heatpng, flank, flanklen):
    ### Int: half of the flank length
    forhalfmer = int(flanklen/2)
    revhalfmer = flanklen - forhalfmer
    forHalfMerSeqs = [''.join(i) for i in itertools.product("ACGT", repeat=forhalfmer)]
    revHalfMerSeqs = [''.join(i) for i in itertools.product("ACGT", repeat=revhalfmer)]
    totalMerSeqs = [''.join(i) for i in itertools.product("ACGT", repeat=flanklen)]
    forMerCount = {}
    revMerCount = {}
    totalMerCount = {}
    for forMerSeq in forHalfMerSeqs:
        forMerCount[forMerSeq] = 0
    for revMerSeq in revHalfMerSeqs:
        revMerCount[revMerSeq] = 0
    for totalMerSeq in totalMerSeqs:
        totalMerCount[totalMerSeq] = 0
    with open(pamcount, "r") as inHandle:
        for line in inHandle:
            if not line.startswith("PAMseq"):
                line = line.strip()
                targetPAMseq, targetPAMcount = line.strip().split("\t")
                targetPAMfor = targetPAMseq[0:forhalfmer]
                targetPAMrev = targetPAMseq[revhalfmer:]
                forMerCount[targetPAMfor] += int(targetPAMcount)
                revMerCount[targetPAMrev] += int(targetPAMcount)
                totalMerCount[targetPAMseq] += int(targetPAMcount)
    
    Generate_heatmap(forMerCount, revMerCount, totalMerCount, heatpng, flank, flanklen)

    return True

def get_file_name_dict(ExpFileList, ConFileList, PAMCountFileTxt):
    resdict = {}
    expcount = 1
    ctrlcount = 1
    with open(PAMCountFileTxt, "w") as outHandle:
        outHandle.write("SampleName\tFileName\tSampleType\n")
        for expfile in ExpFileList:
            sample_name = f"Exp{expcount}"
            resdict[sample_name] = expfile
            outHandle.write(f"{sample_name}\t{os.path.basename(expfile)}\tExp\n")
            expcount += 1
        for ctrlfile in ConFileList:
            sample_name = f"Ctrl{ctrlcount}"
            resdict[sample_name] = ctrlfile
            outHandle.write(f"{sample_name}\t{os.path.basename(ctrlfile)}\tCtrl\n")
            ctrlcount += 1
    return resdict

def Merge_PAM_count(PAMCountList, MergedFile, file_name_dict, upstreamLen):
    ### file_name_dict: {"Exp1": "Exp1.fastq", "Exp2": "Exp2.fastq", "Ctrl1": "Ctrl1.fastq", "Ctrl2": "Ctrl2.fastq"}
    PAMs = [''.join(i) for i in itertools.product("ACGT", repeat=upstreamLen)]
    totalFiles = len(PAMCountList)
    ### Initialize a pandas dataframe
    PAMCountDF = pd.DataFrame(0, index=PAMs, columns=list(file_name_dict.keys()))

    for PAMCountFile in PAMCountList:
        sample_name = None
        base_name = os.path.basename(PAMCountFile).split(".")[0]
        for name, file in file_name_dict.items():
            if base_name in os.path.basename(file):
                sample_name = name
                break
        if sample_name is None:
            raise ValueError(f"File {PAMCountFile} does not match any sample names in the file_name_dict.")
        
        with open(PAMCountFile, "r") as inHandle:
            for line in inHandle:
                if not line.startswith("PAMseq"):
                    line = line.strip()
                    items = line.split("\t")
                    PAMCountDF.at[items[0], sample_name] += int(items[1])
    
    ### Write to the merged file, csv format
    PAMCountDF.to_csv(MergedFile, sep="\t")
    
    return True

def read_and_prepare_data(rawfile):
    data = pd.read_csv(rawfile, header=0, sep="\t")
    data.rename(columns={'Unnamed: 0': 'PAM'}, inplace=True)
    experiment_cols = [col for col in data.columns if col.startswith('Exp')]
    control_cols = [col for col in data.columns if col.startswith('Ctrl')]

    data['LogExperiment'] = np.log2(data[experiment_cols].mean(axis=1) + 1)
    data['LogControl'] = np.log2(data[control_cols].mean(axis=1) + 1)

    data['Difference'] = data['LogExperiment'] - data['LogControl']
    data['DifferenceSmooth'] = data['Difference'].rolling(window=3, min_periods=1).mean()
    
    scaler = RobustScaler()
    data['DifferenceScaled'] = scaler.fit_transform(data['DifferenceSmooth'].values.reshape(-1, 1))
    data['ZscoreScaled'] = zscore(data['DifferenceScaled'])

    _, p_values = ttest_ind(data['LogExperiment'], data['LogControl'], equal_var=False)
    data['p_value'] = p_values


    data['GC_content'] = data['PAM'].apply(lambda x: (x.count('G') + x.count('C')) / len(x))
    data['Length'] = data['PAM'].apply(len)

    return data, experiment_cols, control_cols

def train_model(X_train, y_train):
    clf = RandomForestClassifier()
    param_grid = {
        'n_estimators': [50, 100, 200],
        'max_depth': [None, 10, 20, 30]
    }
    grid_search = GridSearchCV(clf, param_grid, cv=5)
    grid_search.fit(X_train, y_train)
    return grid_search.best_estimator_

def evaluate_model(clf, X_test, y_test, eval_roc_png, eval_prc_png):
    y_pred_proba = clf.predict_proba(X_test)[:, 1]

    fpr, tpr, _ = roc_curve(y_test, y_pred_proba)
    roc_auc = auc(fpr, tpr)
    plt.figure()
    plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic')
    plt.legend(loc="lower right")
    plt.savefig(eval_roc_png)
    plt.close()

    precision, recall, _ = precision_recall_curve(y_test, y_pred_proba)
    plt.figure()
    plt.plot(recall, precision, color='blue', lw=2)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Curve')
    plt.savefig(eval_prc_png)
    plt.close()

def save_classification_report(clf, X_test, y_test, report_file):
    y_pred = clf.predict(X_test)
    report = classification_report(y_test, y_pred, zero_division=0)
    with open(report_file, "w") as f:
        f.write(report)

def generate_logo_data(data, PAMlen):
    significant_data = data[data['ML_significant']]
    logo_data = pd.DataFrame(0.0, index=[str(i+1) for i in range(PAMlen)], columns=['A', 'C', 'G', 'T'])
    totalConfidence = 0
    for index, row in significant_data.iterrows():
        pam = row['PAM']
        if len(pam) != PAMlen:
            print(f"Unexpected PAM length: {pam}")
            continue
        for i, base in enumerate(pam):
            if base in logo_data.columns:
                logo_data.at[f'{i+1}', base] += float(row['Difference'])
                #print(f"Row {index} - Base {base} - Difference {row['Difference']}")
        totalConfidence += row['Confidence']
    

    average_confidence = totalConfidence / len(significant_data)
    overall_confidence_score = (average_confidence - np.min(data['Confidence'])) / (np.max(data['Confidence']) - np.min(data['Confidence']))

    # Scale to [0, 100]
    overall_confidence_score *= 100
    
    for index, row in logo_data.iterrows():
        for base in logo_data.columns:
            if row[base] < 0:   
                logo_data.at[index, base] = abs(row[base])
    
    
    return logo_data, overall_confidence_score

def save_logo_data(logo_data, resfile):
    with open(resfile, "w") as fw:
        fw.write("PO\tA\tC\tG\tT\n")
        for index, row in logo_data.iterrows():
            fw.write(f"{index}\t{row['A']}\t{row['C']}\t{row['G']}\t{row['T']}\n")

def Combine_PAM_score(upstreamPAMscore, downstreamPAMscore, compinedPAMscore, UpstreamLen, DownstreamLen):
    upPAMscore = pd.read_csv(upstreamPAMscore, header=0, sep="\t")
    downPAMscore = pd.read_csv(downstreamPAMscore, header=0, sep="\t")
    totalSize = UpstreamLen + DownstreamLen
    newDataframe = pd.DataFrame(0.0, index=[str(i+1) for i in range(totalSize)], columns=['A', 'C', 'G', 'T'])
    ### Put two df together
    for index, row in upPAMscore.iterrows():
        for base in upPAMscore.columns[1:]:
            newDataframe.at[str(index+1), base] += row[base]
    for index, row in downPAMscore.iterrows():
        for base in downPAMscore.columns[1:]:
            newDataframe.at[str(index+UpstreamLen+1), base] += row[base]

    with open(compinedPAMscore, "w") as fw:
        fw.write("PO\tA\tC\tG\tT\n")
        for index, row in newDataframe.iterrows():
            fw.write(f"{index}\t{row['A']}\t{row['C']}\t{row['G']}\t{row['T']}\n")
    
    return True

def create_pam_logo_combined(compinedPAMscore, combinedPAMlogo, UpstreamLen, DownstreamLen):
    ### UpstreamLen = n, upstreamAnnotation = -n to -1
    ### DownstreamLen = m, downstreamAnnotation = 1 to m

    xannotation = [str(i) for i in range(-int(UpstreamLen), 0)] + [str(i) for i in range(1, int(DownstreamLen)+1)]
    xannotation = ",".join(xannotation)
    logoCMD = "weblogo -f "+compinedPAMscore+" -o "+combinedPAMlogo+" -F jpeg --title combined # 100 --size large --annotate "+xannotation+" --color blue C 'C' --color red T 'T' --color green A 'A' --color orange G 'G' --fineprint CUTPLASMID"
    os.system(logoCMD)

    return True

def create_pam_logo(PAMscoreFile, PAMlogoFile, directions, PAMlen, overall_confidence_score):
    if directions == "upstream":
        xannotation = [str(i) for i in range(-int(PAMlen),-1)]
    else:
        xannotation = [str(i) for i in range(1, int(PAMlen)+1)]
    
    xannotation = ",".join(xannotation)
    logoCMD = "weblogo -f "+PAMscoreFile+" -o "+PAMlogoFile+" -F jpeg --title "+directions+" # "+str(overall_confidence_score)+" --size large --annotate "+xannotation+" --color blue C 'C' --color red T 'T' --color green A 'A' --color orange G 'G' --fineprint CUTPLASMID"
    os.system(logoCMD)

    return True

def visualize_fastp_results(fastp_dir, output_dir):
    fastp_logs = [os.path.join(dp, f) for dp, dn, filenames in os.walk(fastp_dir) for f in filenames if f.endswith(".json")]
    summary_stats = {
        "sample": [],
        "total_reads": [],
        "filtered_reads": [],
        "total_bases": [],
        "filtered_bases": [],
        "q20_rate": [],
        "q30_rate": [],
        "gc_content": []
    }

    for log_file in fastp_logs:
        with open(log_file, "r") as f:
            data = json.load(f)
            sample_name = os.path.basename(log_file).split(".")[0]
            summary_stats["sample"].append(sample_name)
            summary_stats["total_reads"].append(data["summary"]["before_filtering"]["total_reads"])
            summary_stats["filtered_reads"].append(data["summary"]["after_filtering"]["total_reads"])
            summary_stats["total_bases"].append(data["summary"]["before_filtering"]["total_bases"])
            summary_stats["filtered_bases"].append(data["summary"]["after_filtering"]["total_bases"])
            summary_stats["q20_rate"].append(data["summary"]["before_filtering"]["q20_rate"])
            summary_stats["q30_rate"].append(data["summary"]["before_filtering"]["q30_rate"])
            summary_stats["gc_content"].append(data["summary"]["before_filtering"]["gc_content"])

    summary_df = pd.DataFrame(summary_stats)
    summary_df.to_csv(os.path.join(output_dir, "fastp_summary.csv"), index=False)

    plt.figure(figsize=(14, 8))
    plt.subplot(2, 2, 1)
    plt.bar(summary_df["sample"], summary_df["total_reads"], color='blue', alpha=0.7, label='Total Reads')
    plt.bar(summary_df["sample"], summary_df["filtered_reads"], color='red', alpha=0.7, label='Filtered Reads')
    plt.ylabel('Reads')
    plt.title('Total vs Filtered Reads')
    plt.xticks(rotation=90)
    plt.legend()

    plt.subplot(2, 2, 2)
    plt.bar(summary_df["sample"], summary_df["total_bases"], color='blue', alpha=0.7, label='Total Bases')
    plt.bar(summary_df["sample"], summary_df["filtered_bases"], color='red', alpha=0.7, label='Filtered Bases')
    plt.ylabel('Bases')
    plt.title('Total vs Filtered Bases')
    plt.xticks(rotation=90)
    plt.legend()

    plt.subplot(2, 2, 3)
    plt.bar(summary_df["sample"], summary_df["q20_rate"], color='green', alpha=0.7)
    plt.ylabel('Q20 Rate')
    plt.title('Q20 Rate by Sample')
    plt.xticks(rotation=90)

    plt.subplot(2, 2, 4)
    plt.bar(summary_df["sample"], summary_df["q30_rate"], color='orange', alpha=0.7)
    plt.ylabel('Q30 Rate')
    plt.title('Q30 Rate by Sample')
    plt.xticks(rotation=90)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fastp_summary_plots.png"))
    plt.close()

def Process_file(PAMCountFile, PAMscoreFile, PAMlogoFile, evalROC, evalPRC, clfReportFile, PAMlen, directions, tempfile):
    print_status("Reading and preparing data for "+PAMCountFile)
    data, experiment_cols, control_cols = read_and_prepare_data(PAMCountFile)
    if tempfile:
        data.to_csv(tempfile, sep="\t", index=False)

    X = data[['LogExperiment', 'LogControl', 'GC_content', 'Length']]
    
    # Debug: Check the distribution of 'ZscoreScaled' and 'p_value'

    # Adjust condition if necessary
    condition = (data['ZscoreScaled'] < -1.65) & (data['p_value'] < 0.05)
    if not condition.any():
        print("No samples meet the condition, adjusting condition...")
        condition = (data['ZscoreScaled'] < -1) & (data['p_value'] < 0.1)
        if not condition.any():
            print("Still no samples meet the condition, adjusting condition further...")
            condition = (data['ZscoreScaled'] < 0) & (data['p_value'] < 0.5)
    y = condition

    unique_classes = np.unique(y)
    if len(unique_classes) != 2:
        raise ValueError("Expected y to have 2 unique classes, but got: {}".format(unique_classes))

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

    best_clf = train_model(X_train, y_train)
    data['ML_significant'] = best_clf.predict(X)

    y_pred_proba_cv = cross_val_predict(best_clf, X, y, cv=5, method='predict_proba')
    if y_pred_proba_cv.shape[1] == 2:
        data['Confidence'] = y_pred_proba_cv[:, 1]
    elif y_pred_proba_cv.shape[1] == 1:
        data['Confidence'] = y_pred_proba_cv[:, 0]
    else:
        raise ValueError("Unexpected shape of y_pred_proba_cv: {}".format(y_pred_proba_cv.shape))

    evaluate_model(best_clf, X_test, y_test, evalROC, evalPRC)
    save_classification_report(best_clf, X_test, y_test, clfReportFile)

    logo_data, overall_confidence_score = generate_logo_data(data, PAMlen)
    save_logo_data(logo_data, PAMscoreFile)
    create_pam_logo(PAMscoreFile, PAMlogoFile, directions, PAMlen, overall_confidence_score)


### Main function
def main():

    print_header()

    opts = arg_parser()
    ExpFileList = opts.experiment.split(",")
    ExpFileList = [os.path.abspath(x) for x in ExpFileList]
    ConFileList = opts.control.split(",")
    ConFileList = [os.path.abspath(x) for x in ConFileList]
    if len(ExpFileList) == 0 or len(ConFileList) == 0:
        print("Please input the experiment and control data files")
        sys.exit(1)
    OutputDir = opts.output
    OutputDir = os.path.abspath(OutputDir)
    TargetSeq = opts.target
    UpstreamLen = int(opts.upstream)
    DownstreamLen = int(opts.downstream)
    Fastp = opts.fastp
    SeqDirection = opts.seqDirection
    HeatMapFlag = opts.heatmap
    keeptmp = opts.keeptmp

    check_dependencies()
    if not os.path.exists(OutputDir):
        os.makedirs(OutputDir)
    
    if Fastp:
        if not os.path.exists(os.path.join(OutputDir, "fastp/")):
            os.makedirs(os.path.join(OutputDir, "fastp/"))
        if not os.path.exists(os.path.join(OutputDir, "fastp/Exp/")):
            os.makedirs(os.path.join(OutputDir, "fastp/Exp/"))
        if not os.path.exists(os.path.join(OutputDir, "fastp/Ctrl/")):
            os.makedirs(os.path.join(OutputDir, "fastp/Ctrl/"))
        for expFastq in ExpFileList:
            BaseName = os.path.basename(expFastq)
            Fastp_output = os.path.join(OutputDir, "fastp/Exp/", BaseName+"/")
            if not os.path.exists(Fastp_output):
                os.mkdir(Fastp_output)
            Fastp_logfile = os.path.join(Fastp_output, "fastp.log")
            cleanFastp = os.path.join(Fastp_output, BaseName)
            htmlReport = os.path.join(Fastp_output, BaseName+".html")
            jsonReport = os.path.join(Fastp_output, BaseName+".json")
            
            print_status("Filtering the raw data of experiment sample: "+BaseName)
            cmd = "fastp -i "+expFastq+" -o "+cleanFastp+" -h "+htmlReport+" -j "+jsonReport+" > "+Fastp_logfile
            os.system(cmd)
        
        for ctrlFastq in ConFileList:
            BaseName = os.path.basename(ctrlFastq)
            Fastp_output = os.path.join(OutputDir, "fastp/Ctrl/", BaseName+"/")
            if not os.path.exists(Fastp_output):
                os.mkdir(Fastp_output)
            Fastp_logfile = os.path.join(Fastp_output, "fastp.log")
            cleanFastp = os.path.join(Fastp_output, BaseName)
            htmlReport = os.path.join(Fastp_output, BaseName+".html")
            jsonReport = os.path.join(Fastp_output, BaseName+".json")
            
            print_status("Filtering the raw data of control sample: "+BaseName)
            cmd = "fastp -i "+expFastq+" -o "+cleanFastp+" -h "+htmlReport+" -j "+jsonReport+" > "+Fastp_logfile
            os.system(cmd)
        
        ExpFileList = [os.path.join(OutputDir, "fastp/Exp/", os.path.basename(x), x) for x in ExpFileList]
        ConFileList = [os.path.join(OutputDir, "fastp/Ctrl/", os.path.basename(x), x) for x in ConFileList]

        ### Run PAM count
        if not os.path.exists(os.path.join(OutputDir, "PAMcount/")):
            os.makedirs(os.path.join(OutputDir, "PAMcount/"))
        for expFastq in ExpFileList:
            BaseName = os.path.basename(expFastq)
            up_PAM_output = os.path.join(OutputDir, "PAMcount/", BaseName, BaseName+".upstream.PAMcount")
            down_PAM_output = os.path.join(OutputDir, "PAMcount/", BaseName, BaseName+".downstream.PAMcount")
            if not os.path.exists(os.path.join(OutputDir, "PAMcount/", BaseName)):
                os.mkdir(os.path.join(OutputDir, "PAMcount/", BaseName))
            print_status("Counting the PAM of experiment sample: "+BaseName)
            run_PAM_count(expFastq, up_PAM_output, down_PAM_output, TargetSeq, UpstreamLen, DownstreamLen, SeqDirection)
        
        for ctrlFastq in ConFileList:
            BaseName = os.path.basename(ctrlFastq)
            up_PAM_output = os.path.join(OutputDir, "PAMcount/", BaseName, BaseName+".upstream.PAMcount")
            down_PAM_output = os.path.join(OutputDir, "PAMcount/", BaseName, BaseName+".downstream.PAMcount")
            if not os.path.exists(os.path.join(OutputDir, "PAMcount/", BaseName)):
                os.mkdir(os.path.join(OutputDir, "PAMcount/", BaseName))
            print_status("Counting the PAM of control sample: "+BaseName)
            run_PAM_count(ctrlFastq, up_PAM_output, down_PAM_output, TargetSeq, UpstreamLen, DownstreamLen, SeqDirection)
        
        if HeatMapFlag:

            ### Draw k-mer heatmap
            print_status("Drawing k-mer heatmap")
            Heatdir = os.path.join(OutputDir, "PAMheatmap/")
            if not os.path.exists(Heatdir):
                os.makedirs(Heatdir)
            
            for expFastq in ExpFileList:
                BaseName = os.path.basename(expFastq)
                up_PAM_output = os.path.join(OutputDir, "PAMcount/", BaseName, BaseName+".upstream.PAMcount")
                down_PAM_output = os.path.join(OutputDir, "PAMcount/", BaseName, BaseName+".downstream.PAMcount")
                upHeatmap = os.path.join(Heatdir, BaseName+".upstream.heatmap.png")
                downHeatmap = os.path.join(Heatdir, BaseName+".downstream.heatmap.png")
                print_status("Drawing the heatmap of experiment sample: "+BaseName)
                Draw_PAM_Heatmap(up_PAM_output, upHeatmap, "upstream", UpstreamLen)
                Draw_PAM_Heatmap(down_PAM_output, downHeatmap, "downstream", DownstreamLen)
            
            for ctrlFastq in ConFileList:
                BaseName = os.path.basename(ctrlFastq)
                up_PAM_output = os.path.join(OutputDir, "PAMcount/", BaseName, BaseName+".upstream.PAMcount")
                down_PAM_output = os.path.join(OutputDir, "PAMcount/", BaseName, BaseName+".downstream.PAMcount")
                upHeatmap = os.path.join(Heatdir, BaseName+".upstream.heatmap.png")
                downHeatmap = os.path.join(Heatdir, BaseName+".downstream.heatmap.png")
                print_status("Drawing the heatmap of control sample: "+BaseName)
                Draw_PAM_Heatmap(up_PAM_output, upHeatmap, "upstream", UpstreamLen)
                Draw_PAM_Heatmap(down_PAM_output, downHeatmap, "downstream", DownstreamLen)

        ### Merge the PAM count
        UpstreamMergedFile = os.path.join(OutputDir, "PAMcount/UpstreamMergedFile.csv")
        DownstreamMergedFile = os.path.join(OutputDir, "PAMcount/DownstreamMergedFile.csv")
        print_status("Building file-name dictionary for PAM count files")
        PAMCountFileTxt = os.path.join(OutputDir, "PAMcount/PAMCountFile.txt")
        upstreamExpList = [os.path.join(OutputDir, "PAMcount/", os.path.basename(x), os.path.basename(x)+".upstream.PAMcount") for x in ExpFileList]
        downstreamExpList = [os.path.join(OutputDir, "PAMcount/", os.path.basename(x), os.path.basename(x)+".downstream.PAMcount") for x in ExpFileList]
        upstreamCtrlList = [os.path.join(OutputDir, "PAMcount/", os.path.basename(x), os.path.basename(x)+".upstream.PAMcount") for x in ConFileList]
        downstreamCtrlList = [os.path.join(OutputDir, "PAMcount/", os.path.basename(x), os.path.basename(x)+".downstream.PAMcount") for x in ConFileList]
        upstreamTotalList = upstreamExpList + upstreamCtrlList
        downstreamTotalList = downstreamExpList + downstreamCtrlList
        file_name_dict = get_file_name_dict(ExpFileList, ConFileList, PAMCountFileTxt)
        
        print_status("Merge upstream PAM count")
        Merge_PAM_count(upstreamTotalList, UpstreamMergedFile, file_name_dict, UpstreamLen)

        print_status("Merge downstream PAM count")
        Merge_PAM_count(downstreamTotalList, DownstreamMergedFile, file_name_dict, DownstreamLen)

        ### Now define PAM logo storage directory
        LOGOdir = os.path.join(OutputDir, "PAMlogo/")
        if not os.path.exists(LOGOdir):
            os.makedirs(LOGOdir)
        upstreamPAMscore = os.path.join(LOGOdir, "upstream.PAMscore")
        downstreamPAMscore = os.path.join(LOGOdir, "downstream.PAMscore")
        upEvalROC = os.path.join(LOGOdir, "upstream_eval_roc.png")
        upEvalPRC = os.path.join(LOGOdir, "upstream_eval_prc.png")
        downEvalROC = os.path.join(LOGOdir, "downstream_eval_roc.png")
        downEvalPRC = os.path.join(LOGOdir, "downstream_eval_prc.png")
        upPAMLOGO = os.path.join(LOGOdir, "upstream.PAMlogo.jpeg")
        downPAMLOGO = os.path.join(LOGOdir, "downstream.PAMlogo.jpeg")
        upCLFreport = os.path.join(LOGOdir, "upstream.clf_report.txt")
        downCLFreport = os.path.join(LOGOdir, "downstream.clf_report.txt")
        
        ### tempdir
        if keeptmp:
            tempdir = os.path.join(OutputDir, "temp/")
            if not os.path.exists(tempdir):
                os.makedirs(tempdir)

        ### Processing PAM count files
        print_status("Processing upstream PAM count file")
        if keeptmp:
            uptempfile = os.path.join(tempdir, "upstream.PAMcount")
        else:
            uptempfile = None
        Process_file(UpstreamMergedFile, upstreamPAMscore, upPAMLOGO, upEvalROC, upEvalPRC, upCLFreport, UpstreamLen, "upstream", uptempfile)

        if keeptmp:
            downtempfile = os.path.join(tempdir, "downstream.PAMcount")
        else:
            downtempfile = None
        print_status("Processing downstream PAM count file")
        Process_file(DownstreamMergedFile, downstreamPAMscore, downPAMLOGO, downEvalROC, downEvalPRC, downCLFreport, DownstreamLen, "downstream", downtempfile)

        ### Combine PAM score file, for unified visualization
        compinedPAMscore = os.path.join(LOGOdir, "combined.PAMscore")
        combinedPAMlogo = os.path.join(LOGOdir, "combined.PAMlogo.jpeg")

        print_status("Combining PAM score files")
        Combine_PAM_score(upstreamPAMscore, downstreamPAMscore, compinedPAMscore, UpstreamLen, DownstreamLen)

        print_status("Creating combined PAM logo")
        create_pam_logo_combined(compinedPAMscore, combinedPAMlogo, UpstreamLen, DownstreamLen)

        ### Visualize fastp results
        print_status("Visualizing fastp results")
        visualize_fastp_results(os.path.join(OutputDir, "fastp/"), OutputDir)

    else:
        ### Run PAM count
        if not os.path.exists(os.path.join(OutputDir, "PAMcount/")):
            os.makedirs(os.path.join(OutputDir, "PAMcount/"))
        for expFastq in ExpFileList:
            BaseName = os.path.basename(expFastq)
            up_PAM_output = os.path.join(OutputDir, "PAMcount/", BaseName, BaseName+".upstream.PAMcount")
            down_PAM_output = os.path.join(OutputDir, "PAMcount/", BaseName, BaseName+".downstream.PAMcount")
            if not os.path.exists(os.path.join(OutputDir, "PAMcount/", BaseName)):
                os.mkdir(os.path.join(OutputDir, "PAMcount/", BaseName))
            print_status("Counting the PAM of experiment sample: "+BaseName)
            run_PAM_count(expFastq, up_PAM_output, down_PAM_output, TargetSeq, UpstreamLen, DownstreamLen, SeqDirection)
        
        for ctrlFastq in ConFileList:
            BaseName = os.path.basename(ctrlFastq)
            up_PAM_output = os.path.join(OutputDir, "PAMcount/", BaseName, BaseName+".upstream.PAMcount")
            down_PAM_output = os.path.join(OutputDir, "PAMcount/", BaseName, BaseName+".downstream.PAMcount")
            if not os.path.exists(os.path.join(OutputDir, "PAMcount/", BaseName)):
                os.mkdir(os.path.join(OutputDir, "PAMcount/", BaseName))
            print_status("Counting the PAM of control sample: "+BaseName)
            run_PAM_count(ctrlFastq, up_PAM_output, down_PAM_output, TargetSeq, UpstreamLen, DownstreamLen, SeqDirection)
        
        if HeatMapFlag:

            ### Draw k-mer heatmap
            print_status("Drawing k-mer heatmap")
            Heatdir = os.path.join(OutputDir, "PAMheatmap/")
            if not os.path.exists(Heatdir):
                os.makedirs(Heatdir)
            
        for expFastq in ExpFileList:
            BaseName = os.path.basename(expFastq)
            up_PAM_output = os.path.join(OutputDir, "PAMcount/", BaseName, BaseName+".upstream.PAMcount")
            down_PAM_output = os.path.join(OutputDir, "PAMcount/", BaseName, BaseName+".downstream.PAMcount")
            upHeatmap = os.path.join(Heatdir, BaseName+".upstream.heatmap.png")
            downHeatmap = os.path.join(Heatdir, BaseName+".downstream.heatmap.png")
            print_status("Drawing the heatmap of experiment sample: "+BaseName)
            Draw_PAM_Heatmap(up_PAM_output, upHeatmap, "upstream", UpstreamLen)
            Draw_PAM_Heatmap(down_PAM_output, downHeatmap, "downstream", DownstreamLen)
        
        for ctrlFastq in ConFileList:
            BaseName = os.path.basename(ctrlFastq)
            up_PAM_output = os.path.join(OutputDir, "PAMcount/", BaseName, BaseName+".upstream.PAMcount")
            down_PAM_output = os.path.join(OutputDir, "PAMcount/", BaseName, BaseName+".downstream.PAMcount")
            upHeatmap = os.path.join(Heatdir, BaseName+".upstream.heatmap.png")
            downHeatmap = os.path.join(Heatdir, BaseName+".downstream.heatmap.png")
            print_status("Drawing the heatmap of control sample: "+BaseName)
            Draw_PAM_Heatmap(up_PAM_output, upHeatmap, "upstream", UpstreamLen)
            Draw_PAM_Heatmap(down_PAM_output, downHeatmap, "downstream", DownstreamLen)

        ### Merge the PAM count
        UpstreamMergedFile = os.path.join(OutputDir, "PAMcount/UpstreamMergedFile.csv")
        DownstreamMergedFile = os.path.join(OutputDir, "PAMcount/DownstreamMergedFile.csv")
        print_status("Building file-name dictionary for PAM count files")
        PAMCountFileTxt = os.path.join(OutputDir, "PAMcount/PAMCountFile.txt")
        upstreamExpList = [os.path.join(OutputDir, "PAMcount/", os.path.basename(x), os.path.basename(x)+".upstream.PAMcount") for x in ExpFileList]
        downstreamExpList = [os.path.join(OutputDir, "PAMcount/", os.path.basename(x), os.path.basename(x)+".downstream.PAMcount") for x in ExpFileList]
        upstreamCtrlList = [os.path.join(OutputDir, "PAMcount/", os.path.basename(x), os.path.basename(x)+".upstream.PAMcount") for x in ConFileList]
        downstreamCtrlList = [os.path.join(OutputDir, "PAMcount/", os.path.basename(x), os.path.basename(x)+".downstream.PAMcount") for x in ConFileList]
        upstreamTotalList = upstreamExpList + upstreamCtrlList
        downstreamTotalList = downstreamExpList + downstreamCtrlList
        file_name_dict = get_file_name_dict(ExpFileList, ConFileList, PAMCountFileTxt)
        
        print_status("Merge upstream PAM count")
        Merge_PAM_count(upstreamTotalList, UpstreamMergedFile, file_name_dict, UpstreamLen)

        print_status("Merge downstream PAM count")
        Merge_PAM_count(downstreamTotalList, DownstreamMergedFile, file_name_dict, DownstreamLen)

        ### Now define PAM logo storage directory
        LOGOdir = os.path.join(OutputDir, "PAMlogo/")
        if not os.path.exists(LOGOdir):
            os.makedirs(LOGOdir)
        upstreamPAMscore = os.path.join(LOGOdir, "upstream.PAMscore")
        downstreamPAMscore = os.path.join(LOGOdir, "downstream.PAMscore")
        upEvalROC = os.path.join(LOGOdir, "upstream_eval_roc.png")
        upEvalPRC = os.path.join(LOGOdir, "upstream_eval_prc.png")
        downEvalROC = os.path.join(LOGOdir, "downstream_eval_roc.png")
        downEvalPRC = os.path.join(LOGOdir, "downstream_eval_prc.png")
        upPAMLOGO = os.path.join(LOGOdir, "upstream.PAMlogo")
        downPAMLOGO = os.path.join(LOGOdir, "downstream.PAMlogo")
        upCLFreport = os.path.join(LOGOdir, "upstream.clf_report.txt")
        downCLFreport = os.path.join(LOGOdir, "downstream.clf_report.txt")
        
        ### Processing PAM count files
        print_status("Processing upstream PAM count file")
        Process_file(UpstreamMergedFile, upstreamPAMscore, upPAMLOGO, upEvalROC, upEvalPRC, upCLFreport, UpstreamLen, "upstream")

        print_status("Processing downstream PAM count file")
        Process_file(DownstreamMergedFile, downstreamPAMscore, downPAMLOGO, downEvalROC, downEvalPRC, downCLFreport, DownstreamLen, "downstream")

        ### Combine PAM score file, for unified visualization
        compinedPAMscore = os.path.join(LOGOdir, "combined.PAMscore")
        combinedPAMlogo = os.path.join(LOGOdir, "combined.PAMlogo.jpeg")

        print_status("Combining PAM score files")
        Combine_PAM_score(upstreamPAMscore, downstreamPAMscore, compinedPAMscore, UpstreamLen, DownstreamLen)

        print_status("Creating combined PAM logo")
        create_pam_logo_combined(compinedPAMscore, combinedPAMlogo, UpstreamLen, DownstreamLen)

        ### Visualize fastp results
        print_status("Visualizing fastp results")
        visualize_fastp_results(os.path.join(OutputDir, "fastp/"), OutputDir)

    
if __name__ == '__main__':
    main()
