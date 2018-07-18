<!-- README.md is generated from README.Rmd. Please edit that file -->
DMRcompare
==========

This project contains analysis scripts used in the manuscript "An Evaluation of Supervised Methods for Identifying Differentially Methylated Regions in Illumina Methylation Arrays" by Mallik et al. The manual files for all functions in this project is in `docs/DMRcompare.pdf`.

### I. Download Public Methylation Dataset GSE41169

-   `inst/1_downloader2.5.R`

### II. Perform Adjacent Site Clustering to Obtain Clusters of Adjacent CpG Probes (A-clusters)

-   `inst/1_Aclust_data_import.R`

### III. Simulation Study

-   `WriteBumphunterResults()` in script file `R/4_simulate_and_save_Bumphunter_results.R`
-   `WriteDMRcateResults()` in script file `R/4_simulate_and_save_DMRcate_results.R`
-   `WriteProbeLassoResults()` in script file `R/4_simulate_and_save_ProbeLasso_results.R`

The three functions listed above call different functions during each step of the simulation study:

1.  Simulate differentially methylated clusters of CpGs:
    -   `SimulateData()` in script file `R/2_simulatedata.R`
2.  Apply DMR finding methods to the simulated datasets:
    -   `RunBumphunter()` in script file `R/3_RunBumphunter.R`
    -   `RunDMRcate()` in script file `R/3_RunDMRcate.R`
    -   `RunProbeLasso()` in script file `R/3_RunProbeLasso.R`
    -   The `Comb-p` method was implemented in `Python`. The corresponding shell script is `exec/run_combp_working1.sh`
3.  Summarize results of DMR finding methods:
    -   `ProcessBumphunterResults()` in script file `R/5_read_and_summarize_Bumphunter_results.R`
    -   `ProcessDMRcateResults()` in script file `R/5_read_and_summarize_DMRcate_results.R`
    -   `ProcessProbeLassoResults()` in script file `R/5_read_and_summarize_ProbeLasso_results.R`
    -   `ProcessCombpResults()` in script file `R/5_standardize_and_summarize_Comb-p_results.R`

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
