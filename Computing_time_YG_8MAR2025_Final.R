#######################################################################
## Project: Accelerated fitting of joint models with cumulative    ####
##          variations                                             ####
## Script purpose: Numerical examples with cumulative variation    ####
## Date: Mar 8, 2025                                               ####
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

#install.packages("microbenchmark")

#######################################################
## Install R packages                                ##
#######################################################
library(pracma)
library(splines)
library(splines2)
library(microbenchmark)
library(ggplot2)

#######################################################
## Simple numerical examples                         ##
#######################################################

#The time points;
t<-seq(0,10, 0.1)

#S5:
b_1=0.5
b_2=2
b_3=10
b_vec=c(b_1,b_2,b_3)
B_mat<-bs(t,degree=3)
Q_exp<-B_mat%*%b_vec
plot(t, Q_exp,type="l", xlab="t", ylab="Q(t)")

#S6:
b_1=5
b_2=2
b_3=5
b_vec=c(b_1,b_2,b_3)
B_mat<-bs(t,degree=3)
Q_pol<-B_mat%*%b_vec
plot(t,Q_pol,type="l", xlab="t", ylab="Q(t)")

#S7:
b_1=1
b_2=20
b_3=5
b_vec=c(b_1,b_2,b_3)
B_mat<-bs(t,degree=3)
Q_cuv<-B_mat%*%b_vec
plot(t,Q_cuv,type="l", xlab="t", ylab="Q(t)")

#S8:
b_1=10
b_2=2
b_3=2
b_vec=c(b_1,b_2,b_3)
B_mat<-bs(t,degree=3)
Q_for<-B_mat%*%b_vec
plot(t,Q_for,type="l", xlab="t", ylab="Q(t)")

#The plots for all scenarios;
par(mfrow = c(2, 2))
plot(t, Q_exp,type="l", xlab="t", ylab="Q(t)", main="S5")
text(5, 8, expression(Q(t)==0.5*B[1](t)+2*B[2](t)+10*B[3](t) ))
plot(t,Q_pol,type="l", xlab="t", ylab="Q(t)", main="S6")
text(3, 4, expression(Q(t)==5*B[1](t)+2*B[2](t)+5*B[3](t) ))
plot(t,Q_cuv,type="l", xlab="t", ylab="Q(t)",main="S7")
text(6, 2, expression(Q(t)==B[1](t)+20*B[2](t)+5*B[3](t) ))
plot(t,Q_for,type="l", xlab="t", ylab="Q(t)",main="S8")
text(6, 1, expression(Q(t)==10*B[1](t)+2*B[2](t)+2*B[3](t) ))


#######################################################
## Benchmark Analysis for cumulative variation S5-S8 ##
#######################################################

##### Chord;
gt=seq(0,10,by=0.05)
g_grid=bs(gt,degree=3)%*%b_vec 
len_gt<-length(gt)
Chord_fun<-function()
{
  chord<-sum(sqrt((g_grid[-len_gt]-g_grid[-1])^2+0.05^2))
  return(chord)
  
}  
  
  
#### Derivative;  
f1 <- function(s) {sqrt(1+(dbs(s,Boundary.knots=range(t), intercept = FALSE) %*% b_vec)^2)}
Derivative_fun<-function()
  {
  int_res<-integrate(f1, lower=0, upper = 10,abs.tol=1e-03)
  return(int_res)
}


##### Richardson;

f2 <- function(t) c(t, bs(t,degree=3)%*%b_vec)
Richardson_fun<-function()
{
  arc_res<-arclength(f2, 0, 10, tol=1e-03) 
  return(arc_res)
}

set.seed(12345)

mbm <- microbenchmark("Chordal approximation"=Chord_fun(),
                      "Integration + derivative" = Derivative_fun(), 
                      "Richardson integration"=Richardson_fun(),
                      times = 10000)
#View(mbm)
#mbm_S5<-mbm
#write.csv(mbm_S5, "S5_bench.csv")
#save(mbm_S5)

#mbm_S6<-mbm
#write.csv(mbm_S6, "S6_bench.csv")

#mbm_S7<-mbm
#write.csv(mbm_S7, "S7_bench.csv")

#mbm_S8<-mbm
#write.csv(mbm_S8, "S8_bench.csv")


#####################################################
##  Survival Parameters and Data set structure      #
#####################################################
#Constant baseline hazard 
lambda=0.02
B = c(0.05)
X=1
t_val=10

alpha=(0.25)

#######################################################
#######################################################
## Cumulative Hazard time                            ##
#######################################################
#######################################################

#Calculate the cumulative hazard;
#Standard GK nodes and weights;
Kronrod_knots_15= c(-0.949107912342758524526189684047851, -0.741531185599394439863864773280788, 
                    -0.405845151377397166906606412076961, 0,
                    0.405845151377397166906606412076961, 0.741531185599394439863864773280788, 
                    0.949107912342758524526189684047851, -0.991455371120812639206854697526329,
                    -0.864864423359769072789712788640926, -0.586087235467691130294144838258730, 
                    -0.207784955007898467600689403773245, 0.207784955007898467600689403773245,
                    0.586087235467691130294144838258730, 0.864864423359769072789712788640926, 
                    0.991455371120812639206854697526329);
Kronrod_wt_15 = c(0.063092092629978553290700663189204, 0.140653259715525918745189590510238, 
                  0.190350578064785409913256402421014,
                  0.209482141084727828012999174891714, 0.190350578064785409913256402421014, 
                  0.140653259715525918745189590510238, 0.063092092629978553290700663189204,
                  0.022935322010529224963732008058970, 0.104790010322250183839876322541518, 
                  0.169004726639267902826583426598550, 0.204432940075298892414161999234649,
                  0.204432940075298892414161999234649, 0.169004726639267902826583426598550, 
                  0.104790010322250183839876322541518, 0.022935322010529224963732008058970);


#######################################################
## Chordal approximation                             ##
#######################################################
CH_Chord_fun<-function()
{
  
  hz_Arc<-function(t)
  {
  arc_crude<-sum(sqrt((g_grid[-len_gt]-g_grid[-1])^2+0.05^2))
  hz<-lambda*exp(X*B+alpha*arc_crude)
  return(hz)
  }
  
  Hazard = function(t){ 
    wt_scale_vec= (t/2) * Kronrod_wt_15;
    knot_scale_vec = t*(1+Kronrod_knots_15) /2;
    dot(wt_scale_vec, sapply(knot_scale_vec, hz_Arc))
    
  }
  
  #t_val is 10 here;
  Ht<-Hazard(t_val)
  
}  


#######################################################
## Numerical integration method with derivatives     ##
#######################################################

#### Derivative; 
CH_Derivative_fun<-function()
{
        hz_Arc<-function(t) 
        {
         arc_crude<-integrate(f1, lower=0, upper = t,abs.tol=1e-03)
          hz<-lambda*exp(X*B+alpha*as.numeric(arc_crude[1]))
          return(hz)
        }
        
        #hazard rate is a function with respect to t (observed survival time per subject!)
        Hazard = function(t){ 
          #t=10
          wt_scale_vec= (t/2) * Kronrod_wt_15;
          knot_scale_vec = t*(1+Kronrod_knots_15) /2;
          dot(wt_scale_vec, sapply(knot_scale_vec, hz_Arc))
          #hz_Arc(t): this t refer to the knot_scale_vec;
        }
        
        #t_val is 10 here;
        Ht<-Hazard(t_val)

}


#######################################################
## R Arc length with Richardson integration          ##
#######################################################

CH_Richardson_fun<-function()
{
  hz_Arc<-function(t) 
  {
    arc_res<-arclength(f2, 0, t, tol=1e-03) 
    hz<-lambda*exp(X*B+alpha*as.numeric(arc_res[1]))
    return(hz)
  }
  
  #hazard rate is a function with respect to t (observed survival time per subject!)
  Hazard = function(t){ 
    #t=10
    wt_scale_vec= (t/2) * Kronrod_wt_15;
    knot_scale_vec = t*(1+Kronrod_knots_15) /2;
    dot(wt_scale_vec, sapply(knot_scale_vec, hz_Arc))
    #hz_Arc(t): this t refer to the knot_scale_vec;
  }
  
  Ht<-Hazard(t_val)
  
} 

set.seed(12345)
mbm <- microbenchmark("Chordal approximation"=CH_Chord_fun(),
                      "Integration + derivative" = CH_Derivative_fun(), 
                      "Richardson integration"=CH_Richardson_fun(),
                      times = 10000)

par(mfrow = c(2, 2))
autoplot(mbm)

######################################################
## End of Program                                    #
######################################################  


