######################################################################
## Project: Bayesian Lag Joint Model for Survival Outcomes        ####
## Script purpose: Generate the simulated data for the model      ####
## Date: Jul 17st, 2020                                           ####
## Author: Yan Gao (PhD candidate in Biostatistics)               ####
## Division: Epi & Bio, School of Public Health, UIC              ####
## Note: this codes were used to conduct the modeling             ####
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
#install.packages("MASS")
require(MASS)
#install.packages("Matrix")
require(Matrix)
#install.packages("survival")
library(survival)

library(reshape)
detach(package:bayesplot, unload=TRUE)
library(bayesplot)
library(dplyr)
#install.packages("tibble")
library(tibble)
#install.packages("reshape")
library(reshape)

library(reshape2)
#install.packages("stringr")
library(stringr)

library(nlme)
library(lattice)

#install.packages("simstudy")

#library(simstudy)

library(splines)
library(splines2)
library(pracma)

#install.packages("mgcv")
library(mgcv)

#######################################################
## Parallel Processing                               ##
#######################################################
#install.packages("doParallel")
library(doParallel)
n_core<-detectCores()
cl <- parallel::makeCluster(n_core)
doParallel::registerDoParallel(cl)


######################################################
######################################################
## Simulated Bayesian Lag Joint Model Data           #
######################################################
######################################################

#####################################################
## Longitudinal parameter setting                   #
#####################################################
#Total number of subjects
n_sub=10
#totaltimes<-sample(n, min=1, max=10) 
#Fixed total number points for each one
totaltimes<-3
times = seq(0, 12, len = totaltimes)
#Measurement time for the cumulative covariate

#hist(oij)
#View(oij)
#Random intercept and random slopes for each subject
#was formerly the upper triangle
set.seed(123)
n=4
A <- matrix(runif(n^2,0,1)*2, ncol=n) 
Sigma <- t(A) %*% A
c_bs_Sigma<-Sigma
c_bs_Mean<-c(1.2,0.25,0.5,0.9)
#View(c_bs)
#summary(c_bs)
#random error variance
sigma=2

precision<-1/(sigma^2)





#####################################################
##  Longitudinal Data set structure                 #
#####################################################
#Latent variable: random intercept and slope
U<-matrix(,nrow=n_sub,ncol=totaltimes) 
#Observed outcomes with measurement errors
Vij<-matrix(,nrow=n_sub,ncol=totaltimes)


#####################################################
##  Survival Parameters and Data set structure      #
#####################################################
#Constant baseline hazard 
lambda=0.02
B = c(0.05)

alpha=(0.25)
obs.t_mat<-matrix(,n_sub,1)
status_mat<-matrix(,n_sub,1)
id_mat<-matrix(,n_sub,1)

True_t<-matrix(,n_sub,1)
C_t<-matrix(,n_sub,1)

oij=matrix(,n_sub,totaltimes)


#Calculate for the hazard rate, cumulative hazard, survival function;
hazard_dat<-matrix(,n_sub,1)
Cumulative_hazard_dat<-matrix(,n_sub,1)
Survival_dat<-matrix(,n_sub,1)
CV_obs_t_dat<-matrix(,n_sub,1)
#hist(rexp(100,0.05))
#hist(rexp(100,0.3))

#####################################################
##  Data simulation                                 #
#####################################################
#N_sim=1

#Final list to combine all the simulation data!!
Data_final<-list()
        f=1
  
        seed.ID_f<-12345*f
      
        set.seed(seed.ID_f)
        c_bs<-as.matrix(mvrnorm( n_sub,c_bs_Mean, c_bs_Sigma))
        #View(c_bs)
        #summary(c_bs)
        #var(c_bs)
        #Time independent covariate
        set.seed(seed.ID_f)
        X = cbind(age=rnorm(n_sub,40,5)-40)
        #summary(X)
        
        set.seed(seed.ID_f)
        
        #test<-rexp(n_sub*totaltimes, 0.5)
        #summary(test)
        #oij_pre<-matrix(runif(n_sub*totaltimes, 0, 5),byrow=FALSE, ncol=totaltimes)
        #oij_pre<-matrix(0.5+rexp(n_sub*totaltimes, 0.5),byrow=FALSE, ncol=totaltimes)
        oij_pre <- matrix( rep(times, n_sub), byrow=TRUE, ncol=totaltimes)
        
        
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
        
        bound_knots=c(0, 20)
        inner_knots=c(7.5)
        bs_order=3;
        bs_degree=bs_order-1;

        for (i in 1:n_sub) {
          #i=1  
          #ID vector
          print(i)
          id_mat[i,1]=i
          seed.ID_i<-1983*i+107*f

          #####################################################
          ## Survival Data: ti, status, X                     #
          #####################################################
          set.seed(seed.ID_i)
          ui<-runif(1)
   
          c_bs_ID=c_bs[i,]
          X_Cov_ID=X[i]
          
          #Need to return arc_curde, hz, and cumulative hazard and survival functon;
          
          CV_crude<<-function(t) 
          {
            a=0
            n_par=100
            h<-(t-a)/n_par
            x_point<-seq(0,t,by=h)
            total_point<-length(x_point)
            
            #u=5
            
            gu_fun<-function(u)
            {
              dbs_basis<-dbs(u, knots=inner_knots, degree=2, intercept=TRUE, Boundary.knots = bound_knots)[1,]
              c_dot<-c_bs_ID %*% dbs_basis
              gu=sqrt(1+c_dot^2)
              return(gu)
            }
            
            gu_point_all<-sapply(x_point, gu_fun)
            arc_crude<-(2*sum(gu_point_all) - tail(gu_point_all,1) - gu_point_all[1])*h/2
            return(arc_crude)
            
          }
          
          

          hz_Arc<<-function(t) 
          {
                a=0
                n_par=100
                h<-(t-a)/n_par
                x_point<-seq(0,t,by=h)
                total_point<-length(x_point)
                
                #u=5
                
                gu_fun<-function(u)
                {
                    dbs_basis<-dbs(u, knots=inner_knots, degree=2, intercept=TRUE, Boundary.knots = bound_knots)[1,]
                    c_dot<-c_bs_ID %*% dbs_basis
                    gu=sqrt(1+c_dot^2)
                    return(gu)
                }
                
                gu_point_all<-sapply(x_point, gu_fun)
                arc_crude<-(2*sum(gu_point_all) - tail(gu_point_all,1) - gu_point_all[1])*h/2
                hz<-lambda*exp(X_Cov_ID*B+alpha*arc_crude)
            
          }
          
         

          Hazard = function(t){ 
                  wt_scale_vec= (t/2) * Kronrod_wt_15;
                  knot_scale_vec = t*(1+Kronrod_knots_15) /2;
                  CH<-dot(wt_scale_vec, sapply(knot_scale_vec, hz_Arc)  );
                  dot(wt_scale_vec, sapply(knot_scale_vec, hz_Arc)  ) + log(ui) ;
            }
          
          
          Hazard_cal = function(t){ 
            wt_scale_vec= (t/2) * Kronrod_wt_15;
            knot_scale_vec = t*(1+Kronrod_knots_15) /2;
            CH<-dot(wt_scale_vec, sapply(knot_scale_vec, hz_Arc)  );
            return(CH)
           }
          
          
         
          
          
          #Hazard = function(t){ integrate( Vectorize(hz_Arc),0, t, subdivisions=30)$value + log(ui) }
          #integrate function is too slow!!!
          True_ti = uniroot(Hazard,c(0,20),maxiter = 1000)$root
          #View(True_t)
          True_t[i,1]<-True_ti
          
          
          hazard_dat[i,1]<-hz_Arc(True_ti)
          Cumulative_hazard_dat[i,1]<-Hazard_cal(True_ti)
          Survival_dat[i,1]<-exp(-Cumulative_hazard_dat[i,1])
          CV_obs_t_dat[i,1]<-CV_crude(True_ti)
          
          
          #warnings()
          #C_ti<-runif(1,3, Max.FUti)
          set.seed(seed.ID_i)
          C_ti<-2+2*rexp(1,0.02)
          #ti is the observed survival time.
          ti<-pmin(True_ti, C_ti)
          obs.t_mat[i,1]=ti
          #summary(obs.t_mat)
          status<-as.numeric(True_ti<=C_ti)
          status_mat[i,1]=status
     
          #####################################################
          ## Longitudina Data : Vij, oij                      #
          #####################################################
          #The measurement time cannot be larger than the observed survival time!!
          index_i<-(oij_pre[i,]<=obs.t_mat[i,1])
          NA_rep<-rep(NA,length(index_i[index_i == FALSE]))
          sub_timepoint<-length(index_i[index_i == TRUE])
          oij[i,]<-c(oij_pre[i,index_i],NA_rep)

          #oij[i,]<-c(oij_pre[i,])
          for (j in 1:sub_timepoint)
          { 
            #j=2
            #i=500
            U[i,j] <-bs(oij[i,j],knots=inner_knots, degree=2, intercept=TRUE, Boundary.knots = bound_knots)%*%c_bs_ID
            seed.ID<-i*100+10*j+f*12345
            set.seed(seed.ID)
            Vij[i,j]<-rnorm(1, mean=U[i,j], sd=sigma)
          }
          
        }  
        
        #hist(obs.t_mat)

        per_cen<-table(status_mat)[1]/n_sub
        #if (per_cen >0.4) {break}
        #else {
        
        #max(obs.t_mat)
        ######################################################
        ## Survival data: the subject level                 #
        ######################################################
        Dat_surv<-as.data.frame(cbind(id_mat,X,obs.t_mat,status_mat))
        colnames(Dat_surv)<-c("ID","obs_cov","Obs.T","Censor")
        #summary(Dat_surv)
        #table(Dat_surv$Censor)
        #hist(Dat_surv$Obs.T)
        #View(Dat_surv)
        
        ######################################################
        ##Longitudinal data: wide structure subject*timepoint#
        ######################################################
        Dat_long_repeat_time<-as.data.frame(cbind(id_mat,oij))
        #View(Dat_long_repeat_time)
        colnames(Dat_long_repeat_time)<-c("ID","R_T1","R_T2","R_T3")
        #if (miss_1>0.1 | miss_2>0.1) {break}
        #if (miss_e>0.8) {break}
        #hist(Dat_long_repeat_time[,2])
        #hist(Dat_long_repeat_time[,10])
        #Repeated covariates with the cumulative effect
        Dat_long_repeat_cov<-as.data.frame(cbind(id_mat,Vij))
        #save(Dat_long_repeat_cov,file="Long_horizontal_10.RData")
        
        #View(Dat_long_repeat_cov)
        colnames(Dat_long_repeat_cov)<-c("ID","V_T1","V_T2","V_T3")
        
        #summary(Dat_long_repeat_cov)
        #hist(Dat_long_repeat_cov[,2])
        #hist(Dat_long_repeat_cov[,5])
        #hist(Dat_long_repeat_cov[,6])
        #hist(Dat_long_repeat_cov[,10])
        
        ######################################################
        ## Create the longitudinal data                      #
        ######################################################
        cov_long<-melt(Dat_long_repeat_cov, id.vars=c("ID")) %>% arrange(ID) 
        #View(cov_long)
        #names(cov_long)
        cov_long$Time_point<-str_replace(cov_long$variable, "V_", "")
        colnames(cov_long)<-c("ID","Cat","cov","Time_point")
        time_long<-melt(Dat_long_repeat_time, id.vars=c("ID")) %>% arrange(ID)
        time_long$Time_point<-str_replace(time_long$variable, "R_", "")
        colnames(time_long)<-c("ID","Cat2","mt","Time_point")
        #View(time_long)
        cov_long_dat<-inner_join(cov_long,time_long,by=c("ID","Time_point")) %>% 
          arrange(ID, mt) %>%  filter (!is.na(mt))
        
        #View(cov_long_dat)
        
        Data_surv_long<-list(Dat_surv=Dat_surv,cov_long_dat=cov_long_dat, 
                             Dat_long_repeat_time=Dat_long_repeat_time,
                             Dat_long_repeat_cov=Dat_long_repeat_cov,
                             cen_per=per_cen)
        #View(Data_surv_long)
        #Each component of Data_final is a list itself!!
        Data_final[[f]]<-list(Data_surv_long)
        #View(Data_final[[1]])
      
        
        Data_all_comb<-inner_join(Dat_surv, Dat_long_repeat_cov,by=c("ID"))
        
        Data_all_comb_round<-apply(Data_all_comb,2, round,digits=1)
        save(Data_all_comb_round,file="Long_horizontal_10.RData")
        write.csv(Data_all_comb_round, "Data_all_comb_round.csv")
        
        para_all<-cbind(CV_obs_t_dat,hazard_dat,Cumulative_hazard_dat,Survival_dat)
        colnames(para_all)<-c("V_hat","h_hat","H_hat","S_hat")
        para_all_round<-apply(para_all,2, round,digits=1)
        save(para_all_round,file="Para_10.RData")
        write.csv(para_all_round, "para_all_round.csv")
        
        
        
        B_spline_product<-U
        colnames(B_spline_product)<-c("b1B1","b2B2","b3B3")
        save(B_spline_product,file="B_splind_product.RData")
        
        

        #####################################################
        ##  Hypothesis testing for outliers    v            #
        #####################################################
        
        
        
        
        
        
getwd()

######################################################
## End of Program                                    #
######################################################

























