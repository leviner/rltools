# code to read in a single targets export created by Echoview and estimate beam comp 200 kHz channel of 38/200 combi split
# 1) In EV create single target variable using 200 kHz TS. This can be the same as decribed in the Renfree version or using the single-beam method.
# 2) Export the resulting single targets file as CSV and process with this script.
#
# inputs - CSV files (modify file list)
#
# This follows the mthod described in Demer et al., 2015. 
# Use these values to complete an on axis calibration in EV.
# Restrict the single target detection prior to the match ping using the mean of the upper 5% value returned here.
# Then analyze/integrate a 'reduce pings' single target and matched Sv variables. 

library(tidyverse)

rm(list=ls())

# paramters
filenum = 1 # run 1 file at a time, set that here
filterDepth=c(0,6.3,0,0,0) # Set an upper limit line, default should be 0 (There's some shallow crud in the C1_18_19 targets)

max_angle_fit=9 # maximum of-axis angle to allow to fitting routine (e.g. 0.5 Beam width)
parameter_names=c("bw_major", "bw_minor", "offset_major", "offset_minor")
# read in targets
files = c("E:\\MooredEchosounders\\data\\Calibrations\\FinalCalData\\ECSfor200\\C1_17_18_200cal_prefit.target.csv",
            "E:\\MooredEchosounders\\data\\Calibrations\\FinalCalData\\ECSfor200\\C1_18_19_200cal_prefit.target.csv",
            "E:\\MooredEchosounders\\data\\Calibrations\\FinalCalData\\ECSfor200\\C4_17_18_200cal_prefit.target.csv",
            "E:\\MooredEchosounders\\data\\Calibrations\\FinalCalData\\ECSfor200\\C4_18_19_200cal_prefit.target.csv",
            "E:\\MooredEchosounders\\data\\Calibrations\\FinalCalData\\ECSfor200\\C11_17_18_200cal_prefit.target.csv")
file = files[filenum]
targs=read.csv(file)
targs<-subset(targs, Target_range >filterDepth[filenum])
targs_sorted<-sort(targs$TS_uncomp, decreasing=TRUE)
targs_sorted_lin<-10**(targs_sorted/10)
ntargs<-length(targs_sorted_lin)*.05
targs_top5<-targs_sorted_lin[1:ntargs]
x11(width=10, height=10)
par(mfrow=c(2,2))
plot(10*log10(targs_sorted_lin),main=substr(file, 64,71))
abline(h=10*log10(mean(targs_top5)), col="blue")
plot(10*log10(targs_top5),col='red',main=paste("Mean TS of Top 5%:",round(10*log10(mean(targs_top5)),3)))
plot(targs$Target_range, targs$TS_uncomp)
plot(targs$Ping_number, targs$TS_uncomp)
abline(h=10*log10(mean(targs_top5)), col="blue")
print(substr(file, 64,71))
print(10*log10(mean(targs_top5)))
print(10*log10(targs_sorted_lin[1:10]))