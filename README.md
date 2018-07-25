<!-- README.md is generated from README.Rmd. Please edit that file -->
DMRcompare
==========

This project contains analysis scripts used in the manuscript "An Evaluation of Supervised Methods for Identifying Differentially Methylated Regions in Illumina Methylation Arrays" by Mallik et al. The reference files for all functions in this project is in `/docs/DMRcompare.pdf`.

### I. Download Public Methylation Dataset GSE41169

File: `/inst/1_downloader2.5.R`

Input: GEO accession number, criteria for selecting methylation datasets.

Output: `betaVals_mat`, which is beta value matrix for selected methylation samples. This file has rows = cpg ids, columns = sample ids. An example file is at `/data/betaVals_mat.csv`.

### II. Perform Adjacent Site Clustering to Obtain Clusters of Adjacent CpG Probes (A-clusters)

The A-clustering algorithm described in Sofer et al. (2011) was used to identify clusters of adjacent CpGs.

File: `/inst/1_Aclust_data_import.R`

Input:

1.  `betaVals_mat`: a beta value matrix of all CpGs on the array
2.  `cpgLocation_df`: an annotation file that indicates locations of CpGs. This file has rows = cpg ids, columns = chromosome, location. An example file is at `/data/cpgLocation_df.csv`.

Output: `startEndCPG_df`, which is beta value matrix for clusters of CpGs. This file has rows = cpg ids, columns = cluster number, chr, start of cluster, end of cluster, sample ids. An example file is at `/data/startEndCpG_df.csv`.

### III. Simulation Study

There are three main steps in the simulation study. See `/docs/DMRcompare.pdf` for details of each function.

1.  Simulate differentially methylated clusters of CpGs.
    -   File: `SimulateData()` in script file `R/2_simulatedata.R`
    -   Main Input: `betaVals_mat`, `startEndCpG_df` (file that indicates clusters of CpGs), treatment effects to be added to the clusters
    -   Main Output: simulated beta value matrix, where treatment effects were added to 500 randomly selected clusters of CpGs
2.  Apply DMR finding methods to the simulated datasets:
    -   Files:
        -   `RunBumphunter()` in script file `R/3_RunBumphunter.R`
        -   `RunDMRcate()` in script file `R/3_RunDMRcate.R`
        -   `RunProbeLasso()` in script file `R/3_RunProbeLasso.R`
        -   The `Comb-p` method was implemented in `Python`. The corresponding shell script is `exec/run_combp_working1.sh`
    -   Main output: significant DMRs identified by each of the methods. These functions are called by three wrapper functions:
        -   `WriteBumphunterResults()` in script file `R/4_simulate_and_save_Bumphunter_results.R`
        -   `WriteDMRcateResults()` in script file `R/4_simulate_and_save_DMRcate_results.R`
        -   `WriteProbeLassoResults()` in script file `R/4_simulate_and_save_ProbeLasso_results.R`
3.  Summarize results of DMR finding methods:
    -   Files:
        -   `ProcessBumphunterResults()` in script file `R/5_read_and_summarize_Bumphunter_results.R`
        -   `ProcessDMRcateResults()` in script file `R/5_read_and_summarize_DMRcate_results.R`
        -   `ProcessProbeLassoResults()` in script file `R/5_read_and_summarize_ProbeLasso_results.R`
        -   `ProcessCombpResults()` in script file `R/5_standardize_and_summarize_Comb-p_results.R`
    -   Main output: These functions compare the significant DMRs identified by each method, evaluate whether they overlap with the true positive clusters where treatment effects were added, and then compute summary statistics including TP, FP, TN, FN, power, precision, median number of CpGs in significant DMRs

### IV. Table of Results

True Positives (TP), False Positives (FP), False Negatives (FN), Power, Precision, Area under Precision-Recall curve (AuPR), Matthews' correlation coefficient (MCC), F1 Scores (F1) and Elapsed Time (in seconds) for the different DMR detection tools based on simulation datasets:

-   `docs/Method_compare_report_20180705.Rmd`

### V. Plots of Precision and Power over Different Effect Sizes

-   `docs/Method_compare_graphs_20180709.Rmd`

### VI. Plot of Precision and Recall Curves

-   `BuildPRcurve()` in script file `R/6_Build_Precision-Recall_Curve_List.R`
-   `PlotPRCurve()` in script file `R/6_Plot_Precision-Recall_Curves.R`

### VII. Plot of DMR Sizes

-   `docs/Method_compare_graphs_20180709.Rmd`

### VIII. Venn Diagrams for Overlap of DMRs Identified by Each Method

-   `PlotOverlaps()` in script file `R/6_Plot_DMR-Overlaps_Venn.R`
