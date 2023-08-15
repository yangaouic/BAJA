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
require(lattice)
require(mvtnorm)
library(rstan)
library(dplyr )
library(reshape)
library(splines)
library(splines2 )
require(Matrix)
library(lme4)
library(survival)
library(ggplot2)
library(nlme)
library(JM)
#install.packages("rstudioapi")
library(rstudioapi)
#######################################################
## Parallel Processing                               ##
#######################################################
options(mc.cores = parallel::detectCores())
#install.packages("doParallel")
library(doParallel)
n_core<-detectCores()
cl <- parallel::makeCluster(n_core)
doParallel::registerDoParallel(cl)

#######################################################
## Attach the sample data  from JM package           ##
#######################################################
data(aids)
data(aids.id)
pbc2=aids
pbc2.id=aids.id

######################################################
## KM and COX regression - Survival Analysis         #
######################################################

ADS_surv<-pbc2.id
table(ADS_surv$death, ADS_surv$drug)
names(ADS_surv)
surv_dat<-ADS_surv
obs.t<-as.numeric(ADS_surv[,"Time"])
X<-as.numeric(surv_dat[,"drug"])-1 #0: ddC, 1: ddL
event<-surv_dat$death
event_Dat<- surv_dat %>%  filter(death==1)
n.event<-nrow(event_Dat)
obs_t_mat<-as.vector(obs.t)

######################################################
## Data Manipulation for Stan structure             #
######################################################
names(pbc2)
#View(pbc2)
pbc_long<-pbc2 %>% dplyr::select(patient, obstime, CD4) 
pbc_long$cov<-pbc_long$CD4
names(pbc_long)
colnames(pbc_long)<-c("ID", "mt", "CD4","cov")
ADS_long_cov_pre<-pbc_long %>% arrange(ID, mt) %>% dplyr::select(ID, mt, cov)
#View(ADS_long_cov)
#names(ADS_long_cov)

ADS_long_cov<- as.data.frame( apply( ADS_long_cov_pre,2, as.numeric ))
#View(ADS_long_cov)

ADS_long_index_start<-ADS_long_cov %>% mutate(Id_index_start = row_number()) %>%
  group_by(ID) %>% arrange(mt) %>% filter(row_number()==1 ) %>% arrange(ID)
Id_index_start<-ADS_long_index_start$Id_index_start
ADS_long_index_end<-ADS_long_cov %>% mutate(Id_index_end = row_number()) %>%
  group_by(ID) %>% arrange(mt) %>% filter(row_number()==n()) %>% arrange(ID)
Id_index_end<-ADS_long_index_end$Id_index_end
ID_rep_num=Id_index_end-Id_index_start+1


#View(ADS_long_index_start)
#View(ADS_long_index_end)
#View(ADS_long_cov)
Subject<-ADS_long_cov$ID
Measure_time<-ADS_long_cov$mt
#View(Measure_time)
#summary(Measure_time)
long_cov<-ADS_long_cov$cov
n_sub<-length(unique(ADS_long_cov$ID))
n_obs<-length(long_cov)
r_Mat<-cbind(as.vector(rep(1,n_obs)),Measure_time)
Unique_sub_ID<-unique(Subject)


######################################################
## B-spline knots and order and data preparation    ##
## For the calculation of h(t): hazard              ##
## and the likelihood for the longitudinal data     ##
## based on the nonparametric regression model      ##
######################################################

summary(obs.t)
#Calculate the quartiles and choose 0 as the start
#for the B-spline knots.
med_obs_t=median(obs.t)
max_obs_t=max(obs.t)+5
quart_obs_t=c(0, med_obs_t, max_obs_t)
#Pre_specify knots and order for basis functions!!
bs_n_knots=length(quart_obs_t);
bs_order=3;
bs_degree=bs_order-1;
bs_num_basis=bs_n_knots + bs_degree -1; 
#range_t=c(0, End_point)

bound_knots=c(0, max_obs_t)
inner_knots=c(med_obs_t)

bs_inner_knots=inner_knots
range_t=bound_knots

#Calculate the b-spline basis functions results
#for the longitudinal data for all measurement time.
#Measurement time and basis are row and column, respectively
bs_mat_long=bs(Measure_time,knots=bs_inner_knots, 
               degree=bs_degree,intercept =TRUE,
               Boundary.knots = range_t) 
#dim(bs_mat_long) #3949    6
#View(bs_mat_long)

#bs_mat_long<-bs_mat_long[1:7,]
bs_long_block<-t(Matrix::bdiag(split(bs_mat_long,1:nrow(bs_mat_long)))) 
#View(bs_long_block[1:7, 1:20])
#dim(bs_long_block) #  3949 23694=3949*6

sparse_long=extract_sparse_parts(bs_long_block)
long_w=as.vector(sparse_long$w)
long_v=as.vector(sparse_long$v)
long_u=as.vector(sparse_long$u)

long_w_len=length(long_w)
long_v_len=length(long_v)
long_u_len=length(long_u)
######################################################
## Gauss-Kronrod Quadrature knots and weights        #
## For the calculation of H(t): cumulative hazard    #
######################################################

GK_knots=15;
GK_n_obs=n_sub*GK_knots;

#Subject-level
start_t_vec=rep(0, n_sub);
#GK-knopt level
start_knot_vec=rep(0, GK_n_obs);

#Gauss-Kronrod quadrature with Q nodes from R Package JMbayes and Wiki;
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

Kronrod_knots_7=Kronrod_knots_15[1:7];              

Kronrod_wt_7 = c(0.129484966168869693270611432679082, 0.279705391489276667901467771423780, 
                 0.381830050505118944950369775488975, 
                 0.417959183673469387755102040816327, 0.381830050505118944950369775488975, 
                 0.279705391489276667901467771423780, 0.129484966168869693270611432679082);

scal_t=obs_t_mat/2; 

Gauss_choice<-function(num_knots) 
{
  
  if (num_knots==7)
  { knot_vec<<-rep(Kronrod_knots_7,n_sub);
  scal_vec<<-rep(scal_t, each=7);
  wt_vec<<-rep(Kronrod_wt_7, n_sub);
  }
  
  else 
  { knot_vec<<-rep(Kronrod_knots_15,n_sub);
  scal_vec<<-rep(scal_t, each=15);
  wt_vec<<-rep(Kronrod_wt_15, n_sub);
  
  }
  
}


#View(knot_vec)
#View(wt_vec)
#Call the Causs function!!!  
Gauss_choice(GK_knots)

#Calculate the scaled GK knots and weights
#Component-wide product for two vectors: * in R;
#structure: per GK knot per patient for each observation
wt_scale_vec= scal_vec * wt_vec;
knot_scale_vec = rep(obs_t_mat,each=GK_knots)*(1+knot_vec) /2;
#View(wt_scale_vec)
#View(knot_scale_vec)
#range(knot_scale_vec)

######################################################
## Arc Length approximation Setttings                #
## n_seq: number of grids for the arc length         #
######################################################

n_seq=50;  

######################################################
#Subject-level extend to the grid level;     #########
######################################################
Subject_seq=rep(Unique_sub_ID, each=n_seq);

Subject_seq_dat<-as.data.frame(Subject_seq)
#View(Subject_seq_dat)
ADS_seq_index_start<-Subject_seq_dat %>% mutate(Seq_index_start = row_number()) %>%
  group_by(Subject_seq) %>% filter(row_number()==1 ) 
Seq_index_start<-ADS_seq_index_start$Seq_index_start
length(Seq_index_start)
ADS_seq_index_end<-Subject_seq_dat %>% mutate(Seq_index_end = row_number()) %>%
  group_by(Subject_seq) %>%  filter(row_number()==n()) 
Seq_index_end<-ADS_seq_index_end$Seq_index_end
Seq_rep_num=Seq_index_end-Seq_index_start+1


#View(Subject_seq)
mat_rang_t=as.matrix(cbind(start_t_vec,obs_t_mat))
#View(obs_t_mat)
h_t=(obs_t_mat-start_t_vec)/(n_seq-1)
#length(h_t) #500
#View(h_t)
mat_t_seq=t(mapply(seq,mat_rang_t[,1],mat_rang_t[,2],length.out=n_seq))
#View(mat_t_seq)
vec_t_seq=as.vector(t(mat_t_seq))
#View(vec_t_seq)
#per grid point per subject
bs_mat_t=bs(vec_t_seq,knots=bs_inner_knots, 
            degree=bs_degree,intercept =TRUE,
            Boundary.knots = range_t) 
#dim(bs_mat_t) #10000     6
#View(bs_mat_t)

row_bs_t=nrow(bs_mat_t)
col_bs_t=ncol(bs_mat_t)

bs_t_block<-t(Matrix::bdiag(split(bs_mat_t,1:nrow(bs_mat_t)))) 
#dim(bs_t_block) #10000 60000


sparse_t=extract_sparse_parts(bs_t_block)
t_w=as.vector(sparse_t$w)
t_v=as.vector(sparse_t$v)
t_u=as.vector(sparse_t$u)

t_w_len=length(t_w)
t_v_len=length(t_v)
t_u_len=length(t_u)

######################################################
#Subject-level extend to the grid and GK knot level; #
######################################################
Subject_knot=rep(Unique_sub_ID, each=n_seq*GK_knots);
Subject_knot_dat<-as.data.frame(Subject_knot)
#View(Subject_knot_dat)
ADS_knot_index_start<-Subject_knot_dat%>% mutate(Knot_index_start = row_number()) %>%
  group_by(Subject_knot) %>% filter(row_number()==1 ) 

Knot_index_start<-ADS_knot_index_start$Knot_index_start

ADS_knot_index_end<-Subject_knot_dat %>% mutate(Knot_index_end = row_number()) %>%
  group_by(Subject_knot) %>%  filter(row_number()==n()) 

Knot_index_end<-ADS_knot_index_end$Knot_index_end

Knot_rep_num=Knot_index_end-Knot_index_start+1


#View(Subject_knot)
mat_rang_GK=as.matrix(cbind(start_knot_vec,knot_scale_vec))
h_GK=(knot_scale_vec-start_knot_vec)/(n_seq-1)
length(h_GK) #3500
#View(h_GK)
mat_GK_seq=t(mapply(seq,mat_rang_GK[,1],mat_rang_GK[,2],length.out=n_seq))
#View(mat_GK_seq)
vec_GK_seq=as.vector(t(mat_GK_seq))
#per grid point per subject
bs_mat_GK=bs(vec_GK_seq,knots=bs_inner_knots, 
             degree=bs_degree,intercept =TRUE,
             Boundary.knots = range_t) 

#Per grid point per GK knot per subject
#dim(bs_mat_GK) #70000     6


row_bs_GK=nrow(bs_mat_GK)
col_bs_GK=ncol(bs_mat_GK)


bs_GK_block<-t(Matrix::bdiag(split(bs_mat_GK,1:nrow(bs_mat_GK)))) 
dim(bs_GK_block) #70000 420000
sparse_GK=extract_sparse_parts(bs_GK_block)
GK_w=as.vector(sparse_GK$w)
GK_v=as.vector(sparse_GK$v)
GK_u=as.vector(sparse_GK$u)

#View(GK_w)

GK_w_len=length(GK_w)
GK_v_len=length(GK_v)
GK_u_len=length(GK_u)


Tot_grid=n_sub*n_seq
Tot_knot=n_sub*n_seq*GK_knots 

Tot_cs_t=bs_num_basis*Tot_grid
Tot_cs_GK=bs_num_basis*Tot_knot
Tot_cs_all=Tot_cs_t+Tot_cs_GK
######################################################
#In the calculation just call coefficient            #
######################################################

simdata<-list(Subject=Subject,
              Unique_sub_ID=Unique_sub_ID,
              Measure_time=Measure_time,
              r_Mat=r_Mat, 
              long_cov=long_cov,
              n_sub=n_sub,
              n_obs=n_obs,
              obs_t_mat=obs_t_mat,
              X=X,
              event=event,
              bs_n_knots=bs_n_knots,
              bs_order=bs_order,
              bs_num_basis=bs_num_basis,
              quart_obs_t=quart_obs_t,
              bs_mat_long=bs_mat_long,
              
              n_seq=n_seq,
              GK_knots=GK_knots,
              GK_n_obs=GK_n_obs,
              start_t_vec=start_t_vec,
              start_knot_vec=start_knot_vec,
              knot_scale_vec=knot_scale_vec,
              wt_scale_vec=wt_scale_vec,
              Subject_seq=Subject_seq,
              Subject_knot=Subject_knot,
              
              
              #B-splines for the arc length approximation
              h_t=h_t,
              bs_mat_t=bs_mat_t,
              h_GK=h_GK,
              bs_mat_GK=bs_mat_GK,
              
              
              row_bs_t=row_bs_t,
              col_bs_t=col_bs_t,
              
              row_bs_GK=row_bs_GK,
              col_bs_GK=col_bs_GK,
              
              Tot_grid=Tot_grid,
              Tot_knot=Tot_knot, 
              
              Tot_cs_t=Tot_cs_t,
              Tot_cs_GK=Tot_cs_GK,
              Tot_cs_all=Tot_cs_all,
              
              
              #Sparse matrix for the matrix
              long_w=long_w,
              long_v=long_v,
              long_u=long_u,
              
              long_w_len=long_w_len,
              long_v_len=long_v_len,
              long_u_len=long_u_len,
              
              t_w=t_w,
              t_v=t_v,
              t_u=t_u,
              
              t_w_len=t_w_len,
              t_v_len=t_v_len,
              t_u_len=t_u_len,
              
              
              GK_w=GK_w,
              GK_v=GK_v,
              GK_u=GK_u,
              
              GK_w_len=GK_w_len,
              GK_v_len=GK_v_len,
              GK_u_len=GK_u_len,
              
              Id_index_start=Id_index_start,
              Id_index_end=Id_index_end,
              ID_rep_num=ID_rep_num,
              
              Seq_index_start=Seq_index_start,
              Seq_index_end=Seq_index_end,
              Seq_rep_num=Seq_rep_num,
              
              Knot_index_start=Knot_index_start,
              Knot_index_end=Knot_index_end,
              Knot_rep_num=Knot_rep_num
              
)
#search()
#detach("package:Matrix", unload=TRUE)
lowtri=diag(bs_num_basis)
init_fun <- function(...) list(
  sigma=1,
  lambda=1.0,
  beta=1.0, 
  alpha=1.0,
  
  mu=rep(1, bs_num_basis),
  sigma_u=rep(1,bs_num_basis),
  L_u=lowtri, 
  z_u=array(rep(1,bs_num_basis*n_sub),dim=c(bs_num_basis,n_sub)))


#######################################################
## Stan model                                        ##
#######################################################

fitmod <- 
  stan(
    file = "CD4_Yan_Gao_Spline_Stan_2JUN2021.stan",
    data = simdata,
    iter = 5000,
    warmup=1000,
    thin = 1,
    init = init_fun ,
    chains = 4,
    cores = 4,
    control = list(max_treedepth = 10),
    seed=23456)

#View(fitmod)

fit_sum<- as.data.frame(summary(fitmod,pars = c("sigma2","alpha","lambda","beta","cov_mat","mu"))$summary)
fit_sum_c<- as.data.frame(summary(fitmod,pars = c("c_bs"))$summary)

#View(fit_sum_c)
save(fitmod,file="Stan_CD4_Model_CD4.RData")
load("Stan_CD4_Model_CD4.RData")
#View(fit_sum)

PBC_stan_res<-fit_sum[,c("mean","sd","2.5%","50%","97.5%","n_eff","Rhat")]
#View(PBC_stan_res)

PBC_final_dec<- sapply(PBC_stan_res, round, 3)
#View(PBC_final_dec)
write.csv(PBC_final_dec,"CD4_final_Res.csv")

###########################################w##########
##  Extract the posterior data                      ##
######################################################
post_dat<- rstan::extract(fitmod)
#View(post_dat)

###########################################w##########
##  Calculate posterior probability                 ##
######################################################
alpha<-as.data.frame(post_dat$alpha)
#View(alpha)
beta<-as.data.frame(post_dat$beta)
#View(beta)
c_bs.mean <- apply(post_dat$c_bs, c(2,3), mean)

#######################################################
## Figure: Fitting curve for four patients           ##
#######################################################

id_list=c( "87", "390","121","447")
names(pbc2)
#View(pbc2)
pbc_pats_extrem<-pbc2 %>% filter(patient %in% id_list)
pbc_pats_extrem$ID_cha<-paste("Subject", as.character(pbc_pats_extrem$patient))
#View(c_bs.mean)
c_bs.mean_87<-c_bs.mean[,c(87)]
c_bs.mean_390<-c_bs.mean[,c(390)]
c_bs.mean_447<-c_bs.mean[,c(447)]
c_bs.mean_121<-c_bs.mean[,c(121)]

#######################################################
## Chordal approximation                             ##
#######################################################

par(mfrow=c(2,2),mai=c(1, 1, 1,0.8),
    cex.main=1, cex.lab=1, cex.axis=1, cex=1)

#Paitent 87: death: ddC
pat_87<-pbc_pats_extrem %>% filter(patient %in% "87")
t_87<-round(pat_87$Time[1],digits = 2)
t_seq<-seq(0,t_87, by=0.2)
bs_t_seq=bs(t_seq,knots=bs_inner_knots, 
            degree=bs_degree,intercept =TRUE,
            Boundary.knots = range_t) 
pred_87<-bs_t_seq %*% as.vector(c_bs.mean_87)
drug_87=pat_87$drug[1]
plot(pat_87$obstime, pat_87$CD4, ylim=c(0,15), xlim=c(0,20),
     xlab="Time (months)", ylab=expression(sqrt(CD4)), main=paste0("Subject 87"," (",drug_87,")"))
lines(t_seq,pred_87, col="blue",lty = 4,lwd=3)
segments(t_87,-1,t_87,tail(pred_87,1),lwd=2, lty=1, col="red")

gt=seq(0,t_87,by=0.05)
gt_seq=bs(gt,knots=bs_inner_knots, 
          degree=bs_degree,intercept =TRUE,
          Boundary.knots = range_t)
g_grid=gt_seq%*%as.vector(c_bs.mean_87) 
chord<-sum(sqrt((g_grid[-length(gt)]-g_grid[-1])^2+0.05^2))
Gt_mean_87<-round(chord,digits=2)
Gt_mean_87 #10.51
t_87 #10.4

text(10, 13, expression(paste(hat(V)(t), " = 10.51", "; ",  italic(t)," = 10.40",  " (death)"  )))

#Patient 447: ddI, death
pat_447<-pbc_pats_extrem %>% filter(patient %in% "447")
t_447<-round(pat_447$Time[1],digits = 2)
t_seq<-seq(0,t_447, by=0.2)
bs_t_seq=bs(t_seq,knots=bs_inner_knots, 
            degree=bs_degree,intercept =TRUE,
            Boundary.knots = range_t) 
pred_447<-bs_t_seq %*% as.vector(c_bs.mean_447)
drug_447=pat_447$drug[1]
plot(pat_447$obstime, pat_447$CD4, ylim=c(0,15), xlim=c(0,20),
     xlab="Time (months)", ylab=expression(sqrt(CD4)), main=paste0("Subject 447"," (",drug_447,")"))
lines(t_seq,pred_447, col="blue",lty = 4,lwd=3)
segments(t_447,-1,t_447,tail(pred_447,1),lwd=2, lty=1, col="red")
#Cumulative variation estimation by posterior mean of b-spline coefficients!
gt=seq(0,t_447,by=0.05)
gt_seq=bs(gt,knots=bs_inner_knots, 
          degree=bs_degree,intercept =TRUE,
          Boundary.knots = range_t)
g_grid=gt_seq%*%as.vector(c_bs.mean_447) 
chord<-sum(sqrt((g_grid[-length(gt)]-g_grid[-1])^2+0.05^2))
Gt_mean_447<-round(chord,digits=2)
Gt_mean_447 #12.74
t_447 #12.47
text(10, 13, expression(paste(hat(V)(t), " = 12.74", "; ", italic(t), " = 12.47",  " (death)") ))

#Patient 121: ddI, censor
pat_121<-pbc_pats_extrem %>% filter(patient %in% "121")
t_121<-round(pat_121$Time[1],digits = 2)
t_seq<-seq(0,t_121, by=0.2)
bs_t_seq=bs(t_seq,knots=bs_inner_knots, 
            degree=bs_degree,intercept =TRUE,
            Boundary.knots = range_t) 
pred_121<-bs_t_seq %*% as.vector(c_bs.mean_121)
drug_121=pat_121$drug[1]
plot(pat_121$obstime, pat_121$CD4, ylim=c(0,15), xlim=c(0,20),
     xlab="Time (months)", ylab=expression(sqrt(CD4)), main=paste0("Subject 121"," (",drug_121,")"))
lines(t_seq,pred_121, col="blue",lty = 4,lwd=3)
segments(t_121,-1,t_121,tail(pred_121,1),lwd=2, lty=1, col="black")
#Cumulative variation estimation by posterior mean of b-spline coefficients!
gt=seq(0,t_121,by=0.05)
gt_seq=bs(gt,knots=bs_inner_knots, 
          degree=bs_degree,intercept =TRUE,
          Boundary.knots = range_t)
g_grid=gt_seq%*%as.vector(c_bs.mean_121) 
chord<-sum(sqrt((g_grid[-length(gt)]-g_grid[-1])^2+0.05^2))

Gt_mean_121<-round(chord,digits=2)
Gt_mean_121 #12.68 
t_121 #12.50
text(10, 13, expression(paste(hat(V)(t), " = 12.68", "; ", italic(t), " = 12.50",  " (censor)") ))

#Paitent 390: ddC, censor
pat_390<-pbc_pats_extrem %>% filter(patient %in% "390")
t_390<-round(pat_390$Time[1],digits = 2)
t_seq<-seq(0,t_390, by=0.2)
bs_t_seq=bs(t_seq,knots=bs_inner_knots, 
            degree=bs_degree,intercept =TRUE,
            Boundary.knots = range_t) 
pred_390<-bs_t_seq %*% as.vector(c_bs.mean_390)
drug_390=pat_390$drug[1]
plot(pat_390$obstime, pat_390$CD4, ylim=c(0,15), xlim=c(0,20),
     xlab="Time (months)", ylab=expression(sqrt(CD4)), main=paste0("Subject 390"," (",drug_390,")"))
lines(t_seq,pred_390, col="blue",lty = 4,lwd=3)
segments(t_390,-1,t_390,tail(pred_390,1),lwd=2, lty=5, col="black")

gt=seq(0,t_390,by=0.05)
gt_seq=bs(gt,knots=bs_inner_knots, 
          degree=bs_degree,intercept =TRUE,
          Boundary.knots = range_t)
g_grid=gt_seq%*%as.vector(c_bs.mean_390) 
chord<-sum(sqrt((g_grid[-length(gt)]-g_grid[-1])^2+0.05^2))
Gt_mean_390<-round(chord,digits=2)
Gt_mean_390 #12.61
t_390 #12.2
text(6, 8, expression(paste(hat(V)(t), " = 12.61", "; ", italic(t), " = 12.20",  " (censor)") ))

#######################################################
## Pair Plots for Diagnosis                          ##
#######################################################

#pairs(fitmod, pars = c("mu"),las=1)
traceplot(fitmod, pars = c("alpha","lambda","beta","sigma2"), inc_warmup = FALSE)
pairs(fitmod, pars = c("alpha","beta"),las=1)
pairs(fitmod, pars = c("alpha","lambda","beta","sigma2"), las = 1)
pairs(fitmod, pars = c("mu"), las = 1)
pairs(fitmod, pars = c("cov_mat"), las = 1)
pairs(fitmod, pars = c("cov_mat","cor_mat","sigma_u"), las = 1)

pairs(fitmod, pars = c("cor_mat"), las = 1)


pairs(fitmod, pars = c("z_u[1,200]","z_u[2,200]"), las = 1)
traceplot(fitmod, pars = c("mu","sigma2","alpha","lambda","beta"), inc_warmup = FALSE)
traceplot(fitmod, pars = c("sigma_u"), inc_warmup = FALSE)
traceplot(fitmod, pars = c("cov_mat"), inc_warmup = FALSE)
traceplot(fitmod, pars = c("z_u[2,1]"), inc_warmup = FALSE)

plot(fitmod, plotfun = "hist", pars = ("alpha") )
traceplot(fitmod, pars = c("beta"))

traceplot(fitmod, pars = c("sigma_u"))

######################################################
## End of Program                                    #
######################################################  
