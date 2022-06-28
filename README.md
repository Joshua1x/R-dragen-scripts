# R-dragen-scripts
## Scripts in R related to the Illumina DRAGEN-COVID analysis pipeline

### ngs.R
   Main function: select only  patientSamples from the SampleSheet \search for demographic data: lastName, firstName, DOB and sampleCollection date from EpiTrax export (export_15515_XXXXXX.csv)  and/or from the LIMS (NGS_Covid_XXXXX.csv and All_ncovid_NoNGS_to_XXXXX.csv files) and save the ngs file in
/Volumes/NGS/Analysis/covidseq/$run Name


### summary_v2.R
