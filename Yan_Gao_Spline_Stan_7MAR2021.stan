######################################################################
## Project: Bayesian Accelerated Joint Algorithm (BAJA)           ####
## Date: Nov 29th, 2024                                           ####
## Author: Yan Gao, Assistant Professor in Biostatistics          ####
## Division: Division of Biostatistics, MCW                       ####
######################################################################

functions {
  
         //Function to calculate the sum of rows for a matrix;
      row_vector sumRow(matrix inp_dat)
      {
         int num_row=rows(inp_dat);
         row_vector[num_row] y;
      
         for(i in 1:num_row)
         {y[i] = sum(inp_dat[i,]);}
         return y;
      }


      vector Rom_Arc (vector h, int by_Nval,  
                      int Tot_seq,
                      int row_bs,
                      int col_bs,
                      vector bs_basis_w, int[] bs_basis_v, int[] bs_basis_u,
                      vector coeff_vec) 
      {     
      
          int num_t=num_elements(h);
          vector[num_t] h_square;
          matrix[num_t, by_Nval-1] h_sq_vec;

          vector[Tot_seq] g_val_vec;
          matrix[num_t, by_Nval] g_val;
          matrix[num_t,by_Nval-1] g_only;
          matrix[num_t,by_Nval-1] g_duc;
          matrix[num_t,by_Nval-1] g_diff;
          matrix[num_t,by_Nval-1] g_diff_sqare;
          vector[num_t] arc_val;

          h_square= h .* h;
          h_sq_vec=rep_matrix(h_square,by_Nval-1);

          g_val_vec=csr_matrix_times_vector(row_bs, col_bs*row_bs, bs_basis_w, bs_basis_v, bs_basis_u, coeff_vec);
          g_val=(to_matrix(g_val_vec, by_Nval, num_t))';
          g_only=g_val[,1:(by_Nval-1)];
          g_duc=g_val[,2:by_Nval]; 
          g_diff=g_duc-g_only;
          g_diff_sqare=g_diff .* g_diff;
          
          arc_val=to_vector(sumRow(sqrt(h_sq_vec + g_diff_sqare)));
          return(arc_val);
    }

    
    vector exp_cs_fun(matrix rand_ori, int[ ] start_ind, int[ ] end_ind, int[ ] sub_ind_rep, int n_obs)
     {  
        int row_ori=rows(rand_ori);
        int col_ori=cols(rand_ori);
        matrix[row_ori, n_obs] exp_cs; 
        vector[row_ori*n_obs] exp_cs_vec;
        for (n in 1:col_ori) 
          {exp_cs[, start_ind[n]:end_ind[n]] = rep_matrix(rand_ori[,n],sub_ind_rep[n]) ; }
        exp_cs_vec=to_vector(exp_cs);
        return exp_cs_vec;
    } 

}


data {
    int<lower=1> n_obs;            
    real long_cov[n_obs];                
    real<lower=0> Measure_time[n_obs];   
    int<lower=1> n_sub;                   
    int<lower=1,upper=n_sub> Subject[n_obs];  
    vector<lower=0>[n_sub] obs_t_mat; 
    vector[n_sub] X;
    vector<lower=0,upper=1>[n_sub] event;
    vector<lower=1,upper=n_sub>[n_sub] Unique_sub_ID;
    
    int<lower=0> bs_n_knots;
    vector[bs_n_knots] quart_obs_t;
    int<lower=0> bs_num_basis;
    int<lower=0> bs_order;
    matrix[n_obs, bs_num_basis] bs_mat_long;
  
     //Guauss-Kronrod Quadrature Input
     int<lower=1> n_seq;
     int<lower=0> GK_knots;
     int<lower=0> GK_n_obs;
    
     vector<lower=0>[n_sub] start_t_vec;
     vector<lower=0>[GK_n_obs]  start_knot_vec;
     
     vector[GK_n_obs]  knot_scale_vec;
     vector[GK_n_obs]  wt_scale_vec;
     int Subject_seq[n_sub*n_seq];
     int Subject_knot[n_sub*n_seq*GK_knots] ;
     
     vector[n_sub] h_t;
     matrix[n_seq*n_sub, bs_num_basis] bs_mat_t;
     vector[GK_n_obs] h_GK;
     matrix[n_seq*GK_n_obs, bs_num_basis] bs_mat_GK;
     
     int row_bs_t;
     int col_bs_t;
     
     int row_bs_GK;
     int col_bs_GK;
     
     int Tot_grid;
     int Tot_knot;
     
     int Tot_cs_t;
     int Tot_cs_GK;
     int Tot_cs_all;
     
     
     //Sparse matrix associtaed data
     int long_w_len;
     int long_v_len;
     int long_u_len;
     
     vector[long_w_len] long_w;
     int long_v[long_v_len] ;
     int long_u[long_u_len] ;
     
     int t_w_len;
     int t_v_len;
     int t_u_len;
     
     vector[t_w_len] t_w;
     int t_v[t_v_len] ;
     int t_u[t_u_len];
     
     int GK_w_len;
     int GK_v_len;
     int GK_u_len;
     
     vector[GK_w_len] GK_w;
     int GK_v[GK_v_len];
     int GK_u[GK_u_len];
     
     
     int Id_index_start[n_sub];
     int Id_index_end[n_sub];
     int ID_rep_num [n_sub];
     
     int Seq_index_start[n_sub];
     int Seq_index_end[n_sub];
     int Seq_rep_num [n_sub];
     
     int Knot_index_start[n_sub];
     int Knot_index_end[n_sub];
     int Knot_rep_num [n_sub];
    
}


transformed data {
   matrix[n_sub,GK_knots] X_ext_pre;
   vector[GK_n_obs] X_ext;
   X_ext_pre=rep_matrix(X,GK_knots);
   X_ext=to_vector(X_ext_pre');
}

parameters {
  //vector[bs_num_basis] mu;                   
  real<lower=0> sigma;           
  real<lower=0> lambda;
  real beta;
  real alpha;
  
  vector[bs_num_basis] mu;
  
  // random effects standard deviations
  vector<lower=0>[bs_num_basis] sigma_u;       
  // declare L_u to be the Choleski factor of a 2x2 correlation matrix
  cholesky_factor_corr[bs_num_basis] L_u;
  // random effect matrix
  matrix[bs_num_basis,n_sub] z_u;  
  
  //Add the mu parameter in the vector mu[]
  
} 

transformed parameters {
  matrix[bs_num_basis,n_sub]  c_bs;
  matrix[bs_num_basis,n_sub]  mu_exp;
  mu_exp=rep_matrix(mu, n_sub);
  c_bs = mu_exp+diag_pre_multiply(sigma_u, L_u) * z_u;
}

model {
  
  vector[bs_num_basis*n_obs] c_bs_long;
  vector[Tot_cs_all] c_bs_all;
  vector[Tot_cs_t]  c_bs_t; 
  vector[Tot_cs_GK] c_bs_GK; 
  
  vector[n_obs]  u; 
  vector[n_sub]  G_T;
  vector[n_sub]  hT;
  
  vector[GK_n_obs] G_knots_vec;
  vector[GK_n_obs] h_knots_vec;
  
  real  HT;
  
  real  Surv_Log_Lik; 
  
  //c_bs_long is based on the varying number of observations per subject
  c_bs_long=exp_cs_fun(c_bs,Id_index_start,Id_index_end, ID_rep_num,n_obs);
  c_bs_t=exp_cs_fun(c_bs, Seq_index_start, Seq_index_end, Seq_rep_num,  Tot_grid);
  c_bs_GK=exp_cs_fun(c_bs, Knot_index_start, Knot_index_end, Knot_rep_num, Tot_knot);

  u=csr_matrix_times_vector(n_obs, bs_num_basis*n_obs, long_w, long_v, long_u, c_bs_long);
  G_T= Rom_Arc(h_t, n_seq, Tot_grid, row_bs_t, col_bs_t, t_w, t_v, t_u,  c_bs_t);
  hT=lambda*exp(X*beta+alpha*G_T);
            
  G_knots_vec=Rom_Arc(h_GK, n_seq, Tot_knot, row_bs_GK, col_bs_GK, GK_w, GK_v, GK_u, c_bs_GK);
  h_knots_vec=lambda*exp(X_ext*beta+alpha*G_knots_vec);
  //HT: sum over all knots and all subjects and all inside grids;
  HT=dot_product(wt_scale_vec, h_knots_vec);
         
  Surv_Log_Lik = sum(event .* log(hT))-HT;
           
      target += normal_lpdf(long_cov |u, sigma);
      target += Surv_Log_Lik;
      
      //----- Log-priors
      //target += normal_lpdf(sigma_u| 1, 5);
      //target += normal_lpdf(sigma| 1, 5);
      //target += normal_lpdf(mu| 1, 5);
      target += normal_lpdf(lambda| 0.01, 5);
      target += normal_lpdf(beta| 0.01, 10);
      target += normal_lpdf(alpha| 0.01, 10);
      target += lkj_corr_cholesky_lpdf(L_u | 2);
      target += std_normal_lpdf(to_vector(z_u));
     
}


generated quantities {
  real<lower=0> sigma2;
  cov_matrix[bs_num_basis] cor_mat;
  cov_matrix[bs_num_basis] cov_mat;
  
  sigma2=pow(sigma,2);
  // Correlation matrix;
  cor_mat= multiply_lower_tri_self_transpose(L_u); 
  // Variance-Covariance matrix;
  cov_mat= quad_form_diag(cor_mat,sigma_u);
}

/*********************************************************************
**********************************************************************  
*** End of Stan for Simple Case Model on Jan 31, 2021            *****
**********************************************************************  
*********************************************************************/ 



