#!/bin/bash
echo LSB_JOB_REPORT_MAIL=N >> ~/.bashrc
#source /project/bbl_projects/apps/default_modules.sh
#conda activate rstudio
Rscript coupling_accuracy_fx_T.R
