######################################################################
## Project: Bayesian Accelerated Joint Algorithm (BAJA)           ####
## Date: Aug 14th, 2022                                           ####
## Author: Yan Gao, Assistant Professor in Biostatistics          ####
## Division: Division of Biostatistics, MCW                       ####
######################################################################

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
#install.packages("polynom")
library(polynom)
library(splines)
library(splines2)


#Need to use my old PC to run the results!!
par(mfrow=c(1,1))
set.seed(100)
x<-seq(0,10,by=0.05)
length(x)

y_2<-13+x
y_3<-x+x^2
y_4<-5+x+x^2+x^4
y_5<-11*exp(sin(pi*x))

plot(x,y_2)
plot(x,y_3)
plot(x,y_4)
plot(x,y_5)
plot(x, y_2, type="l", col="blue", lwd=3, pch=19, lty=1,
     ylim=c(0,55), xlim=c(0,6), 
     ylab="Q(t)" , xlab="t")
lines(x, y_3, col="red", lty=2,lwd=3)

lines(x, y_4, col="orange", lty=10, lwd=3)
lines(x, y_5, col="black", lty=3, lwd=3)

# Adding a legend inside box at the location (2,40) in graph coordinates.
# Note that the order of plots are maintained in the vectors of attributes.
legend("topright",legend=c("S1: Q(t) = 13 + t",expression("S2: Q(t) = "~ t+t^2),
                          expression("S3: Q(t) = 5+ "~ t+t^2+t^4),
                     "S4: Q(t) = 11*exp{sin(pi*t)}"), 
       col=c("blue","red","orange","black"),
       lty=c(1,2,10,3), ncol=1,  bty = "n")


#Arc length
f1<-function(t) {c(t,13+t)}
arclength(f1,a=0, b=10) #14.14214

f2<-function(t) {c(t,t+t^2)}
arclength(f2,a=0, b=10) # 110.7356

f3<-function(t) {c(t,5+t+t^2+t^4)}
arclength(f3,a=0, b=10) # 10110.26

f4<-function(t) {c(t,11*exp(sin(pi*t)))}
arclength(f4,a=0, b=10) #259.1568

#arclength: tol = 1e-05, nmax=20

arclength


#######################################################
## Convergence: S1 and S2 with closed-form solution  ##
#######################################################

gf1<-function(t) {13+t}
gf2<-function(t) {t+t^2}
gf3<-function(t) {5+t+t^2+t^4}
gf4<-function(t) {11*exp(sin(pi*t))}

#######################################################
## Chordal approximation convergence plots           ##
#######################################################



#######################################################
## S1                                                ##
#######################################################

num_grid_all<-seq(10, 1000, by=10)
grid_count<-length(num_grid_all)

#### For gf1: gf1 is a straight line thus no estimation error.
chord_mat<-matrix(,nrow=grid_count, ncol=2)
for (i in 1:grid_count) {
  #i=1
  num_grid<-num_grid_all[i]
  gt=seq(0,10,length.out=num_grid)
  g_grid=gf1(gt) 
  delta=10/(num_grid-1)
  chord<-sum(sqrt((g_grid[-length(gt)]-g_grid[-1])^2+delta^2))
  res_chord=c(num_grid,chord)
  chord_mat[i,]=res_chord
  
}

chord_mat_gf1<-chord_mat
#Converge immediately because S1 is a line;
par(mfrow=c(2,2))
plot(chord_mat[,1], chord_mat[,2], xlab="Number of grids", ylab="Estimated Cumulative Variation",
     ylim=c(14.12, 14.16), main="Convergence plot for S1", xlim=c(10,1000))
abline(h=sqrt(2)*10, col="blue", lwd=3, lty=1)

#true value: sqrt(2)*10



#######################################################
## S2                                                ##
#######################################################

##### For gf2: gf2 is a curve
chord_mat<-matrix(,nrow=grid_count, ncol=2)

for (i in 1:grid_count) {
  #i=1
  num_grid<-num_grid_all[i]
  gt=seq(0,10,length.out=num_grid)
  g_grid=gf2(gt) 
  delta=10/(num_grid-1)
  chord<-sum(sqrt((g_grid[-length(gt)]-g_grid[-1])^2+delta^2))
  res_chord=c(num_grid,chord)
  chord_mat[i,]=res_chord
  
}


#110.7356 at number of grides=250;
chord_mat_gf2<-chord_mat

plot(chord_mat[,1], chord_mat[,2], xlab="Number of grids", ylab="Estimated Cumulative Variation",
     ylim=c(110.7, 110.75),  main="Convergence plot for S2")
#axis(1, at = seq(10, 1000, by=50), las=2)

t_val_1<-10
t_val_2<-1+2*t_val_1
t_val_3<-0
t_val_4<-1+2*t_val_3

gf2_true<-( (log (sqrt(1+t_val_2^2) + t_val_2) + t_val_2*sqrt(1+t_val_2^2)) -
           (log (sqrt(1+t_val_4^2) + t_val_4) + t_val_4*sqrt(1+t_val_4^2)) )/4


abline(h=gf2_true, col="blue", lwd=3, lty=1)



#######################################################
## S3                                                ##
#######################################################

##### For gf3: No closed form
chord_mat<-matrix(,nrow=grid_count, ncol=2)

for (i in 1:grid_count) {
  #i=1
  num_grid<-num_grid_all[i]
  gt=seq(0,10,length.out=num_grid)
  g_grid=gf3(gt) 
  delta=10/(num_grid-1)
  chord<-sum(sqrt((g_grid[-length(gt)]-g_grid[-1])^2+delta^2))
  res_chord=c(num_grid,chord)
  chord_mat[i,]=res_chord
  
}


#Grid 80 converge to 10110.26;
chord_mat_gf3<-chord_mat

plot(chord_mat[,1], chord_mat[,2], xlab="Number of grids", ylab="Estimated Cumulative Variation", 
     main="Convergence plot for S3", ylim=c(10110.16,10110.28 ))

#######################################################
## S4                                                ##
#######################################################

##### For gf4: No closed form
chord_mat<-matrix(,nrow=grid_count, ncol=2)

for (i in 1:grid_count) {
  #i=1
  num_grid<-num_grid_all[i]
  gt=seq(0,10,length.out=num_grid)
  g_grid=gf4(gt) 
  delta=10/(num_grid-1)
  chord<-sum(sqrt((g_grid[-length(gt)]-g_grid[-1])^2+delta^2))
  res_chord=c(num_grid,chord)
  chord_mat[i,]=res_chord
  
}


chord_mat_gf4<-chord_mat
#converge to  259.14 at grid=810;

plot(chord_mat[,1], chord_mat[,2], xlab="Number of grids",  ylab="Estimated Cumulative Variation",
     main="Convergence plot for S4", ylim=c(135,265)  )




######################################################
## End of Program                                    #
######################################################  



