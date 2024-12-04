######################################################################
## Project: Bayesian Accelerated Joint Algorithm (BAJA)           ####
## Date: Nov 29th, 2024                                           ####
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
library(rstan)
library(rstan)
library(dplyr)
library(reshape)
library(splines)
library(splines2)
require(Matrix)
require(lattice)
require(mvtnorm)
#library(JMbayes)


#######################################################
## Attach the sample data PBC2 from JMbayes          ##
#######################################################

load("pbc2.RData")
load("pbc2.id.RData")

names(pbc2)

pbc2$status2 <- as.numeric(pbc2$status != "alive")
pbc2.id$status2 <- as.numeric(pbc2.id$status != "alive")

#View(pbc2)

pbc_pats_extrem<-pbc2 %>% filter(id %in% id_list)
pbc_pats_extrem$ID_cha<-paste("Subject", as.character(pbc_pats_extrem$id))

#View(pbc_pats_extrem)
my.settings <- canonical.theme(color=FALSE)
my.settings[['strip.background']]$col <- "gray"
my.settings[['strip.border']]$col<- "white"

lattice::xyplot(log(serBilir) ~ year |ID_cha, data = pbc_pats_extrem,
                #groups=ID,
                subset = id %in% id_list, type = c("p"), 
                lwd = 2, layout = c(3, 2),
                xlab=list(label="Time (years)", fontsize=12),
                par.settings=my.settings ,
                ylab=list(label="log(serum Bilirubin)", fontsize=14),
                panel=function(...) {
                  panel.xyplot(...)
                  panel.grid()
                }
)


######################################################
## KM and COX regression - Survival Analysis         #
######################################################

ADS_surv<-pbc2.id
names(ADS_surv)
surv_dat<-ADS_surv
obs.t<-as.numeric(ADS_surv[,"years"])
X<-as.numeric(surv_dat[,"age"])
event<-surv_dat$status2
event_Dat<- surv_dat %>%  filter(status2==1)
n.event<-nrow(event_Dat)
obs_t_mat<-as.vector(obs.t)

#View(ADS_surv)
#summary(ADS_surv)

######################################################
## Data Manipulation for Stan structure             #
######################################################
names(pbc2)
#View(pbc2)
pbc_long<-pbc2 %>% select(id, year,serBilir) 
pbc_long$cov<-log(pbc_long$serBilir)
names(pbc_long)
colnames(pbc_long)<-c("ID", "mt", "SerBilir","cov")
ADS_long_cov_pre<-pbc_long %>% arrange(ID, mt) %>% select(ID, mt, cov)
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

n_seq=20;  

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

#options(mc.cores = parallel::detectCores()) 
fitmod <- 
  stan(
    file = "Yan_Gao_Spline_Stan_7MAR2021.stan",
    data = simdata,
    iter = 2000,
    thin = 1,
    init = init_fun ,
    chains = 1,
    cores = 1,
    #refresh=100,
    #control = list(adapt_delta = 0.99),
    control = list(max_treedepth = 15),
    seed=12345)

fit_sum<- as.data.frame(summary(fitmod,pars = c("sigma2","alpha","lambda","beta","cov_mat","mu"))$summary)


fit_sum_c<- as.data.frame(summary(fitmod,pars = c("c_bs"))$summary)

View(fit_sum_c)
save(fitmod,file="Stan_PBC_Model_2.RData")
load("Stan_PBC_Model_2.RData")

View(fit_sum)

PBC_stan_res<-fit_sum[,c("mean","sd","2.5%","50%","97.5%","n_eff")]

View(PBC_sta_res)

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
#View(fitmod)
#print(fitmod, pars = c("mu"), digits = 3)



#######################################################
## Fitting Curves                                    ##
#######################################################


#id_list=c("56","187","93", "30","217", "20","144","2")

id_list=c("56","187","93", "30","217", "20","144","2", "134", "38","51","70")



pbc_pats_extrem<-pbc2 %>% filter(id %in% id_list)
pbc_pats_extrem$ID_cha<-paste("Subject", as.character(pbc_pats_extrem$id))

names(pbc_pats_extrem)
#View(pbc_pats_extrem)

View(pbc2)


pars <- rstan::extract(fitmod)
#View(pars)
c_bs<-pars$c_bs[1,,]
c_bs.mean <- apply(pars$c_bs, c(2,3), mean)

#View(c_bs.mean)



c_bs.mean_56<-c_bs.mean[,c(56)]
c_bs.mean_187<-c_bs.mean[,c(187)]
c_bs.mean_93<-c_bs.mean[,c(93)]
c_bs.mean_30<-c_bs.mean[,c(30)]
c_bs.mean_217<-c_bs.mean[,c(217)]
c_bs.mean_20<-c_bs.mean[,c(20)]

c_bs.mean_2<-c_bs.mean[,c(2)]
c_bs.mean_144<-c_bs.mean[,c(144)]


c_bs.mean_51<-c_bs.mean[,c(51)]
c_bs.mean_70<-c_bs.mean[,c(70)]

t_seq<-seq(0,19, by=0.2)
bs_t_seq=bs(t_seq,knots=bs_inner_knots, 
               degree=bs_degree,intercept =TRUE,
               Boundary.knots = range_t) 


#Plot B-spline basis functions!
par(mfrow=c(2,2))
plot(bs_t_seq[,1]~t_seq, ylim=c(0,max(bs_t_seq)), type='l', lwd=2, col=1, 
     xlab="B-spline basis", ylab="B1")

plot(bs_t_seq[,2]~t_seq, ylim=c(0,max(bs_t_seq)), type='l', lwd=2, col=2, 
     xlab="B-spline basis", ylab="B2")

plot(bs_t_seq[,3]~t_seq, ylim=c(0,max(bs_t_seq)), type='l', lwd=2, col=3, 
     xlab="B-spline basis", ylab="B3")


plot(bs_t_seq[,4]~t_seq, ylim=c(0,max(bs_t_seq)), type='l', lwd=2, col=4, 
     xlab="B-spline basis", ylab="B4")


#for (j in 2:ncol(bs_t_seq)) lines(bs_t_seq[,j]~t_seq, lwd=2, col=j)

#legend("topleft", legend=c("B1", "B2", "B3", "B4"),
       #col=c("black","red", "green","blue"), lty=1:4,cex=0.8)

#col 1: black
#2:red
#3: green
#4: blue




#library(ggplot2)




pbc_pats_extrem<-pbc2 %>% filter(id %in% id_list)
pbc_pats_extrem$ID_cha<-paste("Subject", as.character(pbc_pats_extrem$id))



my.settings <- canonical.theme(color=FALSE)
my.settings[['strip.background']]$col <- "gray"
my.settings[['strip.border']]$col<- "white"

lattice::xyplot(log(serBilir) ~ year |ID_cha, data = pbc_pats_extrem,
                #groups=ID,
                subset = id %in% id_list, type = c("p"), 
                lwd = 2, layout = c(3, 2),
                xlab=list(label="Time (years)", fontsize=12),
                par.settings=my.settings ,
                ylab=list(label="log(serum Bilirubin)", fontsize=14),
                panel=function(...) {
                  panel.xyplot(...)
                  panel.grid()
                }
)





par(mfrow=c(2,5))


pat_2<-pbc_pats_extrem %>% filter(id %in% "2")
pat_2$log_ser<-log(pat_2$serBilir)
pred_2<-bs_t_seq %*% as.vector(c_bs.mean_2)
plot(pat_2$year, pat_2$log_ser, ylim=c(-1,5), xlim=c(0,20),
     xlab="Time (year)", ylab="log(serum Bilirubin)", main="Subject 2")
lines(t_seq,pred_2, col="blue",lty = 1)



pat_20<-pbc_pats_extrem %>% filter(id %in% "20")
pat_20$log_ser<-log(pat_20$serBilir)
pred_20<-bs_t_seq %*% as.vector(c_bs.mean_20)
plot(pat_20$year, pat_20$log_ser, ylim=c(-1,5), xlim=c(0,20),
     xlab="Time (year)", ylab="log(serum Bilirubin)", main="Subject 20")
lines(t_seq,pred_20, col="blue",lty = 1)

pat_30<-pbc_pats_extrem %>% filter(id %in% "30")
pat_30$log_ser<-log(pat_30$serBilir)
pred_30<-bs_t_seq %*% as.vector(c_bs.mean_30)
plot(pat_30$year, pat_30$log_ser, ylim=c(-1,5), xlim=c(0,20),
     xlab="Time (year)", ylab="log(serum Bilirubin)", main="Subject 30")
lines(t_seq,pred_30, col="blue",lty = 1)

pat_51<-pbc_pats_extrem %>% filter(id %in% "51")
pat_51$log_ser<-log(pat_51$serBilir)
pred_51<-bs_t_seq %*% as.vector(c_bs.mean_51)
plot(pat_51$year, pat_51$log_ser, ylim=c(-1,5), xlim=c(0,20),
     xlab="Time (year)", ylab="log(serum Bilirubin)", main="Subject 51")
lines(t_seq,pred_51, col="blue",lty = 1)


pat_56<-pbc_pats_extrem %>% filter(id %in% "56")
pat_56$log_ser<-log(pat_56$serBilir)

pred_56<-bs_t_seq %*% as.vector(c_bs.mean_56)
plot(pat_56$year, pat_56$log_ser, ylim=c(-1,5), xlim=c(0,20),
     xlab="Time (year)", ylab="log(serum Bilirubin)", main="Subject 56")
lines(t_seq,pred_56, col="blue",lty = 1)



pat_70<-pbc_pats_extrem %>% filter(id %in% "70")
pat_70$log_ser<-log(pat_70$serBilir)
pred_70<-bs_t_seq %*% as.vector(c_bs.mean_70)
plot(pat_70$year, pat_70$log_ser, ylim=c(-1,5), xlim=c(0,20),
     xlab="Time (year)", ylab="log(serum Bilirubin)", main="Subject 70")
lines(t_seq,pred_70, col="blue",lty = 1)


pat_93<-pbc_pats_extrem %>% filter(id %in% "93")
pat_93$log_ser<-log(pat_93$serBilir)
pred_93<-bs_t_seq %*% as.vector(c_bs.mean_93)
plot(pat_93$year, pat_93$log_ser, ylim=c(-1,5), xlim=c(0,20),
     xlab="Time (year)", ylab="log(serum Bilirubin)", main="Subject 93")
lines(t_seq,pred_93, col="blue",lty = 1)



pat_144<-pbc_pats_extrem %>% filter(id %in% "144")
pat_144$log_ser<-log(pat_144$serBilir)
pred_144<-bs_t_seq %*% as.vector(c_bs.mean_144)
plot(pat_144$year, pat_144$log_ser, ylim=c(-1,5), xlim=c(0,20),
     xlab="Time (year)", ylab="log(serum Bilirubin)", main="Subject 144")
lines(t_seq,pred_144, col="blue",lty = 1)


pat_187<-pbc_pats_extrem %>% filter(id %in% "187")
pat_187$log_ser<-log(pat_187$serBilir)
pred_187<-bs_t_seq %*% as.vector(c_bs.mean_187)
plot(pat_187$year, pat_187$log_ser, ylim=c(-1,5), xlim=c(0,20),
     xlab="Time (year)", ylab="log(serum Bilirubin)", main="Subject 187")
lines(t_seq,pred_187, col="blue",lty = 1)



pat_217<-pbc_pats_extrem %>% filter(id %in% "217")
pat_217$log_ser<-log(pat_217$serBilir)
pred_217<-bs_t_seq %*% as.vector(c_bs.mean_217)
plot(pat_217$year, pat_217$log_ser, ylim=c(-1,5), xlim=c(0,20),
     xlab="Time (year)", ylab="log(serum Bilirubin)", main="Subject 217")
lines(t_seq,pred_217, col="blue",lty = 1)






######################################################
## End of Program                                    #
######################################################  

