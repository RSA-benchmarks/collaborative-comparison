#!R
#=======================================================
# author: J.A.Postma
# contact: digigram at gmail dot com
# date: summer 2010
# license: gpl
# purpose: Simulating nutrient uptake by roots
# reference: Itoh S, Barber SA. 1983. A numerical solution of whole plant nutrient uptake for soil-root systems with root hairs. Plant and Soil 70: 403-413.
# content: example code for using the model. 
# requirements: the barber-cushman-onlineVersion0.1-function-library.r in the working dir of R
# disclaimer: This software comes without any warrenty, and may contain errors.
#
#=======================================================

#make sure your working directory is set correctly, use getwd() and setwd() to change it. 
#load model functions
source("barber-cushman-onlineVersion0.32-function-library.r")

#set parameter set
#param=roothairPublication1983
param=paramBarberBookP

#param is a list
str(param)

#change and param with
param$b=163

#simple run, first argument the parmeter set to be used, second a simple label for output message (mainly useful when running batch jobs)
m<-barber(param,msg="my test run")

#plot the results
plotbarber(param,m)

#integrate the results, in case you want numbers instead of graphs
res<-barberTotN(param,m)
names(res)

#do a sensitivity analysis (may take a long time, depending on the list or parameters!)
dat<-sensitivityBarber(param, sensitivity_names=c("Imax","Km","De"))
#("Cli","r0","r1","Imax","Km","De","b","v0","Imaxh","lh","rh","Nh","Nh~lh"))
plotSensitivityBarber(dat)

#write sensitivity analysis to file
plotSensitivityBarber(dat,filename="sensitivityAnalysis")

#test functions
const<-setConstants(param) #check that the constants in test are right
plotCn(param) #check that the differentiation functions are right



