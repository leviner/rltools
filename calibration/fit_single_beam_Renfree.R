# Code to read in a single targets export created by Echoview and estimate beam compensation on 200 kHz channel of 38/200 combi split
# 1) In EV create single target variable with 38  kHz angles and 200 kHz TS
# 2) Set the ecs file to have a 359.9 deg beam angle for major and minor beam axes for the 200 kHz (to effectively turn off beam comp during single targe detection)
# 3) Export the resulting single targets file as CSV and process with this script.
#
# inputs - CSV file
# outputs -  beam_fit_results - fitted [bw_major, bw_minor, offset_major, offset_minor]
#
# This the same as the method described in  Renfree et al, 2019
# but allows the beam to be fit in both along and athwarship, and I keep the coordiante system in the 38 kHz's frame of reference
# The approach assumes that 38 kHz and 200 kHz are at the same range as max separation is say 20 cm ((6^2+.25^2)^.5)=6.005, so small angle approx applies.
#
# Originally produced by A. De Robertis
# Modified by R. Levine

library(tidyverse)

rm(list=ls())

# paramters
max_angle_fit=9 # maximum of-axis angle to allow to fitting routine (e.g. 0.5 Beam width)
parameter_names=c("bw_major", "bw_minor", "offset_major", "offset_minor")
# read in targets
targs=read.csv("E:\\MooredEchosounders\\data\\Calibrations\\FinalCalData\\ECSfor200\\C1_18_19_200cal_prefit.target.csv") # 1049

# set paramters [bw_major, bw_minor, offset_major, offset_minor]
par_angles=c(18, 18, 0, 0)

#look at headers
glimpse(targs)
targs=mutate(targs,Off_axis_angle=(Angle_major_axis^2+Angle_minor_axis^2)^.5)

#  a quick plot of the targtets
ggplot(data=targs,mapping=aes(x=Angle_major_axis,y=Angle_minor_axis,color=TS_uncomp))+geom_point()+scale_color_gradientn(colours = rainbow(10))

# now filter targets to a narrower distributions
targs_cropped<<-filter(targs,Off_axis_angle<max_angle_fit) # assingn as a global variable (the <<-)
ggplot(data=targs_cropped,mapping=aes(x=Angle_major_axis,y=Angle_minor_axis,color=TS_uncomp))+geom_point()+scale_color_gradientn(colours = rainbow(10))

# add a function to do the beam compensation
TS_fit <- function(par_angles, targs_cropped) {
  x=(2*targs_cropped$Angle_major_axis-par_angles[3])/par_angles[1]
  y=(2*targs_cropped$Angle_minor_axis-par_angles[4])/par_angles[2]
  comp=6.0206*(x^2+y^2 -(.18*x^2*y^2))
  TS=targs_cropped$TS_uncomp+comp
  # now calculate the mean squared error in compensation
  # root mean squared errror from average of all measurements
  # error measures mean square difference in TS across all measurements from mean value 
  # i.e. if TS is not dependent on angle, then error will be low
  error=(mean((TS-mean(TS))^2))^.5
  return(error) }


#a = TS_fit(par_angles = par_angles, targs_cropped = targs_cropped)
beam_fit_results=optim(par = par_angles, fn = TS_fit, targs_cropped = targs_cropped)
parameter_names
beam_fit_results[[1]]
rms_error=TS_fit(beam_fit_results[[1]],targs_cropped)
rms_error