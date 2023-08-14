#######################################################################
## Project: Accelerated fitting of joint models with cumulative    ####
##          variations                                             ####
## Script purpose: Numerical examples with cumulative variation    ####
## Date: Aug 3, 2023                                               ####
## Author: Yan Gao, Assistant Professor in Biostatistics           ####
## Division: Division of Biostatistics, MCW                        ####
## Note: this codes were used to illustrate computational gains    ####
#######################################################################

######################################################
####      Clean the temporary files          #########
######################################################
#Clean the console
cat("\014") 
dir()
#clean all objects
rm(list = ls())
#System information
sessionInfo()
getwd()

#######################################################
## Install R packages                                ##
#######################################################
library(pracma)
library(splines)
library(splines2)

#######################################################
## Simple numerical examples                         ##
#######################################################

#The time points;
t<-seq(0,10, 0.1)

#Scenario 1:
b_1=0.5
b_2=2
b_3=10
b_vec=c(b_1,b_2,b_3)
B_mat<-bs(t,degree=3)
Q_exp<-B_mat%*%b_vec
plot(t, Q_exp,type="l", xlab="t", ylab="Q(t)")

#Scenario 2:
b_1=5
b_2=2
b_3=5
b_vec=c(b_1,b_2,b_3)
B_mat<-bs(t,degree=3)
Q_pol<-B_mat%*%b_vec
plot(t,Q_pol,type="l", xlab="t", ylab="Q(t)")

#Scenario 3:
b_1=1
b_2=20
b_3=5
b_vec=c(b_1,b_2,b_3)
B_mat<-bs(t,degree=3)
Q_cuv<-B_mat%*%b_vec
plot(t,Q_cuv,type="l", xlab="t", ylab="Q(t)")

#Scenario 4:
b_1=10
b_2=2
b_3=2
b_vec=c(b_1,b_2,b_3)
B_mat<-bs(t,degree=3)
Q_for<-B_mat%*%b_vec
plot(t,Q_for,type="l", xlab="t", ylab="Q(t)")

#The plots for all scenarios;
par(mfrow = c(2, 2))
plot(t, Q_exp,type="l", xlab="t", ylab="Q(t)", main="Scenario 1")
text(5, 8, expression(Q(t)==0.5*B[1](t)+2*B[2](t)+10*B[3](t) ))
plot(t,Q_pol,type="l", xlab="t", ylab="Q(t)", main="Scenario 2")
text(3, 4, expression(Q(t)==5*B[1](t)+2*B[2](t)+5*B[3](t) ))
plot(t,Q_cuv,type="l", xlab="t", ylab="Q(t)",main="Scenario 3")
text(6, 2, expression(Q(t)==B[1](t)+20*B[2](t)+5*B[3](t) ))
plot(t,Q_for,type="l", xlab="t", ylab="Q(t)",main="Scenario 4")
text(6, 1, expression(Q(t)==10*B[1](t)+2*B[2](t)+2*B[3](t) ))

#######################################################
## Numerical integration method with derivatives     ##
#######################################################

start = Sys.time()
for (i in 1:10000) {
  f <- function(s) { sqrt(1+(dbs(s,Boundary.knots=range(t), intercept = FALSE) %*% b_vec)^2)}
  int_res<-integrate(f, lower=0, upper = 10,abs.tol=1e-03)
} 
time_int<-Sys.time() - start
int_ave<-time_int/10000
int_ave

#Scenario 1:
#value: 14.89 , time:0.0004 secs
#Scenario 2:
#value: 11.58 , time:0.0004 secs
#Scenario 3:
#value: 165.47 , time:0.0030 secs
#Scenario 4:
#value: 13.67 , time:0.0012 secs


#######################################################
## R Arc length with Richardson integration          ##
#######################################################

start = Sys.time()
for (i in 1:10000) {
f <- function(t) c(t, bs(t,degree=3)%*%b_vec)
arc_res<-arclength(f, 0, 10, tol=1e-03) 
} 
time_arc<-Sys.time() - start
arc_ave<-time_arc/10000
arc_ave

#Scenario 1:
#value: 14.89 , time:0.0013 secs
#Scenario 2:
#value: 11.58, time:0.0009 secs
#Scenario 3:
#value: 165.47 , time:0.0016 secs
#Scenario 4:
#value: 13.67 , time: 0.0010 secs


#######################################################
## Chordal approximation                             ##
#######################################################
start = Sys.time()
for (i in 1:10000) {
gt=seq(0,10,by=0.05)
g_grid=bs(gt,degree=3)%*%b_vec 
chord<-sum(sqrt((g_grid[-length(gt)]-g_grid[-1])^2+0.05^2))
} 
our_arc<-Sys.time() - start
our_ave<-our_arc/10000
our_ave
#Scenario 1:
#value: 14.89, time:0.0002 secs
#Scenario 2:
#value: 11.58, 0.0002 secs, 
#Scenario 3:
#value: 165.47, time:0.0002  secs
#Scenario 4:
#value: 13.67 , time: 0.0002 secs


cho<-c(0.0002,0.0002,0.0002,0.0002)
Int<-c(0.0004,0.0004,0.0030,0.0012)
Ric<-c(0.0013,0.0009,0.0016,0.0010)

cho/Int
cho/Ric

######################################################
## End of Program                                    #
######################################################  


