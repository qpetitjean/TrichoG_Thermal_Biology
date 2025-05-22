# TrichoG_Thermal_Biology
## General
![GitHub](https://img.shields.io/github/license/qpetitjean/TrichoG_Thermal_Biology)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/qpetitjean/TrichoG_Thermal_Biology)
![GitHub last commit](https://img.shields.io/github/last-commit/qpetitjean/TrichoG_Thermal_Biology)
![GitHub repo size](https://img.shields.io/github/repo-size/qpetitjean/TrichoG_Thermal_Biology)

Datasets: TODO


Permanent copy of this code repository: TODO


This repository contains:
  * Some R useful R function to allow some computations, and reproduce the visualization of the results within the **R_Func** directory.
  * The scripts used to process tracking data (clean, append temperature data, compute ans smooth movement metrics, display change in movement metrics along temperature changes) within the **TrajectoriesProcessing** directory.
  * The scripts to test Trichogramma strains' sensitivity to temperature within the **Analyses** directory.
  * Image processing macro (imageJ) to detect arena contours and create distance maps within the **OtherTools** directory.

Files within this code repository, Datas from Zenodo repository and the step by step procedure detailed below allow to fully reproduce the results of the paper .....

These results are discussed in the following manuscript:


## Step by step procedure to reproduce the analyses:

### Raw video tracking processing (Optional)

As storing video files and raw video-tracking output is ressource consuming we have uploaded the raw video-tracking output for only 1 assay (both warm and cold ramp for `2021-08-06-AM002`). However, the following step by step procedure has been used to analyse all video-recorded assays. This section is thus mostly made to depict the raw results processing pipeline. Another, although shorter example (truncated warm ramp), can be find within the `How to` section of the MoveR R package website at the following URL:
https://qpetitjean.github.io/MoveR/articles/MoveR-Clean-FilterData.html and https://qpetitjean.github.io/MoveR/articles/MoveR-ComputeMetrics.html


So, let's start to process the raw data (video-tracking output). To do it, run the scripts in the following order: 

1- CleaningScript.r # clean the raw tracking output by: 
                                  * filtering out Infinite values (lost individuals)
                                  * filtering out detected particles with length outside the 95% CI of the individuals length                
                                  * compute speed of individual and filter out values above the 999th quantile
                                  * compute distance to the edge and filter out trajectory part detected outside (in case more than half of the mean body length of individual is outside the arena)

The script automatically save the results in a **Results** directory. 
It includes:
  * an histogram of the frequency of infinite value detected before filtering (`Infvalues.svg`)
  * an histogram of log10 transformed individual length and 95% CI (`IndLen.svg`)
  * an histogram of log10 transformed individual speed (bimodal distribution) and 999th quantile (`IndSpeed.svg`)
  * an histogram of the distance to the edge, with positive values corresponding to individuals outside the arena (`IndOut.svg`)
  * an histogram of log10 transformed trajectory length after filtration steps (`TrackLen.svg`)
  * a summary of the amount of data (number of trajectories, length of the trajectories, percent of data kept) kept after each filtration steps (`FilterSummary.csv`)
  * summaries of the video and trajectories characteristics (video length, number and length of trajectories) before (`VideoSummaryBeforeFilter`) and after filtration steps (`VideoSummaryAfterFilter`).
  * the cleaned dataset for further use (`cleaned_[AssayID].csv`)


2- TempAppendScript.r # append the measured temperature data to the cleaned tracking data
The script automatically save the results in a **Results** directory. 
It includes:
  * a plot displaying the evolution the measured temperature along time (`Temp_Ramp.svg`)
  * the cleaned dataset with temperature ramp appended to each trajectory (`Clean_Temp_[AssayID].csv`)


3- MetricsComputScript.r # compute and smooth various metrics (e.g., speed, activity, sinuosity) on the cleaned tracking data using parallel processing
The script automatically save the results in a **Results** directory. 
It includes:
  * an histogram of log10 transformed speed and automatically detected threshold classifying behavioral state as active or inactive (`SpeedactivTresh.svg`)
  * two log files returning infos about parallel processing (`logMetricCompute.txt` and `logMetricSmooth.txt`)
  * the cleaned dataset with computed/smoothed metrics (`Clean_Temp_metrics_[AssayID].csv`)


4- ActivityCompScript.r # identify activity states (active vs. inactive) using a 2d array with turning angle variance and log10 transformed speed over trajectories
The script automatically save the results in a **Results** directory. 
It includes: 
  * an Interactive 3D density plot of activity states (`3d_activityClustPlot.html`)
  * the Warm and cold datasets with activity classification appended (`Clean_Temp_metrics_activ_[AssayID].csv`)
  

5- SmoothMetricScriptBIDIME.r # computes smoothed movement metrics over time using a 90s sliding window sampled every 5s, weighted by tracklet length
The script automatically saves the results in the Results directory.
It includes:
  * a log file returning infos about parallel processing (`logSmooth.txt`)
  *  the raw output of the smoothing step (`temporalTrend` R function from `MoveR`) including dataset of all metric appended and smoothed over the timeline (`Clean_Temp_metrics_activ_smoothed_[AssayID]_outputWtdALL.rds`)
  * the smoothed metrics and temperature data (`Clean_Temp_metrics_activ_smoothed_[AssayID].csv`)
  * the smoothed temporal trends of each metrics (15 .svg file named `SmoothedMetrics_[metric].svg`)
  

6- CIBootstrapScript.r # computes 95% bootstrapped (student) confidence intervals (CI) around smoothed movement metrics using a 90s sliding window sampled every 5s
The script automatically saves the results in the Results directory.
It includes:
  * a log file returning infos about parallel processing (`logSmooth.txt`)
	* the raw output of the bootstrapping step (`temporalBoot` R function from `MoveR`), including Studentized confidence intervals for each metric (`Clean_Temp_metrics_activ_smoothed95CI_[AssayID]_outputWtdBootsALL.rds`)
	* a final CSV dataset containing 95% CI and mean values per time point for all metrics (`Clean_Temp_metrics_activ_smoothed95CI_[AssayID].csv`)
	* 15 .svg files displaying temporal trends with 95% confidence envelopes per metric (`SmoothedMetrics95CI_[metric].svg`)






7- CorrectBootstrapOutliersScriptBIDIME.r # optional, allows to remove spurious CI estimation (performed on a low number of trajectories/fragment or outside a specified -expected- range)
8- FinalPlotScriptBIDIME.r # ensure that you have the Strains_plotv3.r file, allow to plot the computed metric, it is better to not plot more than 3 metrics on the same plot
9- FinalMetricScriptBIDIME.r # allow to extract various predefined metrics such as global mean activity, CTMAX, CTMIN....

Other: 
To plot repeated ramp and check recovery by visual inspection run:
RepetPlot.r script # ensure that you also have the Strains_plotv4-Repetability.r file

To plot both warm and cold ramps in a two-panel plot run:
FinalPlotPairedScriptBIDIME.r # ensure that you have the Strains_plotv3_pairedRamp.r file


### Perform statistical analyses on extracted metrics

