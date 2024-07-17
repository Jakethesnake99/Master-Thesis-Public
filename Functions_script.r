
# 0. libraries needed (might be double loaded but who cares)
library(tidyverse) # used for data wrangling
library(sampling) # used for sampling
library(MASS) # used for generating data
library(e1071) # used for CV 
library(caret) # used for CV 
library(parallel) # used for making clusters on cores
library(doParallel) # used for parallelizing the simulation 


# 1. Functions used in the analysis
## 1.1 This function generates a population of size N with p covariates with an even covariance of "multi_cor" will break if multi_cor is too high!!!!!
population_func = function(N,p,multi_cor){
  
  # make a covariance matrix for p variables + 1 (selectivity variable) with a variance of 1
  cov_matrix = diag(rep(0, p + 1))
  cov_matrix[1:(p), 1:(p)] = multi_cor
  cov_matrix[p+1, 1:(p+1)]  = 0
  cov_matrix[1:(p+1),(p+1)]  = 0
  diag(cov_matrix) = 1
  
  #create p +1 variables with a 0 mean and above cov_matrix
  Pop_df = as.data.frame(mvrnorm(n = N, mu = rep(0,p+1), Sigma = cov_matrix))
  
  #Rename variables and shuffle df
  Pop_df  = Pop_df %>%
    rename_with(~ paste0("X", seq_along(.)), 1:p) %>%
    rename(z = p + 1) %>%
    sample_n(nrow(.)) %>%
    return(Pop_df)
}

## 1.2 Takes as coefficient vector (betas) of length p and generates a new vector (phis) of length p which when added to the coefficient vector (betas) results in a set correlation and bias factor relative to betas (w)
adj_bias_func = function(betas, target_correlation, bias_factor, p) {
  set.seed(2135248)
  
  #sets max size (bias) of new vector 
  max_sum = sum(betas) * abs(bias_factor)
  
  #Takes a vector and rescales it and checks its correlation
  obj_func = function(new_vec) {
    phi = new_vec / sum(new_vec) * max_sum
    
    if (bias_factor > 0) {
      w =  betas + phi  
    } else {
      w =betas - phi  
    }

    corr = cor(betas, w)
    return(abs(corr - target_correlation))
  }
  
  #generates a ranom vector as start point and then uses optim to rotate it to desired correlation 
  initial_vec = runif(p, 0.00, max_sum)
  result = optim(par = initial_vec, fn = obj_func, method = "L-BFGS-B", lower = 0.00, upper = max_sum)
  phi = result$par / sum(result$par) * max_sum
  return(phi)
}


## 1.3 This functions generates a zero-centered outcome. The outcome is based on the above population function with p X variables and one additional z variable which is interacted with all X variables. Betas is a combined vector of true betas + the additional vector phi. y_var sets the residual variance
y_func = function(Pop_df,betas,y_var){
  
  #checks how many X variables have been generated 
  num_cols = ncol(Pop_df)

  for (i in 1:(num_cols - 1)) {
    # Create interaction variable between the current column and the last one
    interaction_col = paste0("X", i, "_interaction")
    Pop_df[[interaction_col]] = Pop_df[,i] * Pop_df[,num_cols]
  }
  
  
  #Generate errors and a model matrix with all relevant variables 
  error=rnorm(nrow(Pop_df),0,y_var)
  pop_matrix = model.matrix(~. -1 -z, data = Pop_df)
  
  Pop_df = Pop_df %>%
    mutate(y = pop_matrix %*% betas + error) %>%
    relocate(y)
  
  return(Pop_df)
}


## 1.4 transform and rescales the selecitivity column function
Selectivity_func = function(column,nps_n,scaling = 2,shifting = 1){
  column_out = 1/(1 + exp(scaling*-(column-shifting)))
  column_out = nps_n*column_out  / sum(column_out)
}


## 1.5 Baysian support functions (SEE WISLONWSKY paper for explanation of code)
calc.posterior = function(mu_0, k_0, Vin = diag(length(mu_0)), a_0 = 0, b_0 = 0,
                          y = Y, x = X) {
  if (dim(as.matrix(x))[2] !=  length(mu_0) | length(y) !=  dim(as.matrix(x))[1])
    print("Dimensions mismatch!") else {
      n_mu_0 = length(mu_0)
      p = dim(as.matrix(x))[2]
      n = length(y)
      V = Vin * k_0
      # some OLS values
      mu_hat = as.matrix(lm(y ~ x - 1)$coefficients)
      xtx = t(x) %*% x
      RSS = t(y - x %*% mu_hat) %*% (y - x %*% mu_hat)
      # posterior mean of mu
      Sigma_t = solve(xtx + solve(V))
      W = Sigma_t %*% xtx
      mu_mean = W %*% mu_hat + (diag(p) - W) %*% as.matrix(mu_0)
      SS = RSS + t(mu_hat - as.matrix(mu_0)) %*% solve(solve(xtx) + V) %*%
        (mu_hat - as.matrix(mu_0))
      v = n + 2 * a_0
      mu_var = Sigma_t * (as.numeric(SS) + 2 * b_0)/(n + 2 * a_0 - 2)
      # posterior of tau
      tau_mean = (n/2 + a_0)/(SS/2 + b_0)
      tau_var = (n/2 + a_0)/(SS/2 + b_0)^2
      return(list(mu_mean = mu_mean, mu_cov = sqrt(diag(mu_var)), tau_mean_sd = c(tau_mean,
                                                                                  sqrt(tau_var))))
    }
}

fun.hot.c = function(hot.n,n){if (hot.n < 0.05) 1/log(n) else 1/n}

hotelling.test  = function(lm_prob, lm_np){
  xbar = lm_prob$coefficients
  mu_0 = lm_np$coefficients
  vcovm = vcov(lm_prob)
  p = length(lm_prob$coefficients)
  n = length(lm_prob$residuals)
  t2 = t(xbar - mu_0)%*%solve(vcovm)%*%(xbar - mu_0)
  f = p*(n-1)/(n-p)
  Fstat = t2/f
  p_val = pf(Fstat,p,n-p,lower.tail = F)
  return(p_val)
}

gen.res = function(PS_y,ps_mat, nps_y, nps_mat){
  #sample_prob - sample of prob data
  #sample of nprob data
  # conjugate non-inf result
  yp_std = PS_y
  xp_std = ps_mat
  coefsize = dim(xp_std)[2] #number of coefficients
  k = length(yp_std) # length of Prob sample
  
  result_ni = calc.posterior(mu_0 = rep(0,coefsize), k_0 = k, y = yp_std, x = xp_std)
  #ML prob data
  lm_pobj = lm(yp_std ~ xp_std-1)
  # nonprob data
  ynp_std = nps_y
  xnp_std = nps_mat
  nnp = length(ynp_std) #np sample length
  #ML for nonprob
  lm_npobj = lm(ynp_std ~ xnp_std-1)
  #hotelling test
  hot.n = hotelling.test(lm_pobj,lm_npobj)
  #conjugate posterior #k_0fun.hot.c(hot.n,nnp)
  result_c = calc.posterior(mu_0 = lm_npobj$coefficients,
                            k_0 = fun.hot.c(hot.n,nnp), y = yp_std,x = xp_std)
  
  result = as.array(result_c$mu_mean)
  return(result)
}


## 1.6 Cross-validation for optimal hyper parameter function
CV_func = function(PS_DF, w, folds) {
  
  # set n and the range of lambda values to search through 
  n = nrow(PS_DF)
  lambda = c(seq(0.0005, 0.99, length.out = 150), seq(1,5,length.out = 100),seq(6, 100, length.out = 95),10^5)
  
  #do cv to find best lambda for the RRE, DBTLE and ABTLE 
  CV_results_inv = rowMeans(sapply(folds, function(fold) {
    train_df = PS_DF[-fold, ]
    test_df = PS_DF[fold, ]
    test_y = test_df$y
    train_y = train_df$y
    train_mat = model.matrix(y ~ .-1, data = train_df)
    test_mat = model.matrix(y ~ .-1, data = test_df)
    
    I  = diag(ncol(train_mat))
    A = t(train_mat) %*% train_mat
    B = t(test_mat) %*% test_y
    C = t(test_mat) %*% test_mat
    D =  t(train_mat) %*% train_y
    
    sapply(lambda, function(lam) {
      
      J = solve(A + (n * lam *I))
      above = t(J %*% w) %*% B - t(J %*% D) %*% C %*% J %*% w
      denom = n * t(J %*%w) %*% C %*% (J %*% w)
      
      eta = as.numeric(above/denom)
      if(eta<0){
        eta=0
      }
      # Compute beta_hat_* for all ridge variants (ABTLE, RRE, DBTLE)
      beta_hat_AB = J %*% (D + n * eta * w)
      beta_hat_RRE = J %*% D
      beta_hat_DB = J %*% (D + n * lam * w)
      # Calculate mean squared error
      mse_RRE=mean((test_mat %*% beta_hat_RRE - test_y)^2)
      mse_AB=mean((test_mat %*% beta_hat_AB - test_y)^2)
      mse_DB=mean((test_mat %*% beta_hat_DB - test_y)^2)
      
      return(c(mse_AB,mse_RRE,mse_DB))
    })
  }))
  #find best value and extract the corresponding lambda value
  AB_vec = CV_results_inv[seq(1, length(CV_results_inv), by = 3)]
  RRE_vec = CV_results_inv[seq(2, length(CV_results_inv), by = 3)]
  DB_vec = CV_results_inv[seq(3, length(CV_results_inv), by = 3)]
  
  min_ind = which(AB_vec ==  min(AB_vec))
  lambda_optimal_AB = lambda[min_ind]
  min_ind = which(RRE_vec ==  min(RRE_vec))
  lambda_optimal_RRE = lambda[min_ind]
  
  min_ind = which(DB_vec ==  min(DB_vec))
  lambda_optimal_DB = lambda[min_ind]
  
  #only for the ABTLE repeat the step but using only the previously found lambda value 
  results = sapply(folds, function(fold) {
    train_df = PS_DF[-fold, ]
    test_df = PS_DF[fold, ]
    test_y = test_df$y
    train_y = train_df$y
    train_mat = model.matrix(y ~ .-1, data = train_df)
    test_mat = model.matrix(y ~ .-1, data = test_df)
    I  = diag(ncol(train_mat))
    
    
    A = t(train_mat) %*% train_mat
    B = t(test_mat) %*% test_y
    C = t(test_mat) %*% test_mat
    D =  t(train_mat) %*% train_y
    
    J = solve(A + (n * lambda_optimal_AB *I))
    above = t(J %*% w) %*% B - t(J %*% D) %*% C %*% J %*% w
    denom = n * t(J %*%w) %*% C %*% (J %*% w)
    
    
    return(list(above = above, denom = denom))
  })
  # rather than looking at each fold sum them up to get \eta^hat
  total_above  = sum(unlist(results[1,]))
  total_denom  = sum(unlist(results[2,]))
  
  # Divide the total above by the total denom
  eta_optimal_AB = total_above / total_denom
  if(eta_optimal_AB<0){
    eta_optimal_AB=0
  }
  #return all penalty parameters
  return(list(lambda_optimal_AB = lambda_optimal_AB, eta_optimal_AB = eta_optimal_AB,lambda_optimal_RRE=lambda_optimal_RRE,lambda_optimal_DB=lambda_optimal_DB))
}


## 1.6 analysis function
analysis_func = function(Pop_df,sample_col,ps_n, nps_n,main_coef,r,popres23) {
  
  # sample PS and NPS and select the subset of relevant variables
  Pop_df = Pop_df %>%
    mutate(
      PS = ifelse(row.names(.) %in% sample(nrow(.), ps_n), 1, 0),
      NPS = UPrandomsystematic(Pop_df[[sample_col]])
    )
  NPS_DF = subset(Pop_df, NPS ==  1)[,1:(p+1)]
  PS_DF = subset(Pop_df, PS ==  1)[,1:(p+1)]
  
  
  #get coefficients for each variable and the perform CV (above function))
  w_hat = coef(lm(y ~ .-1, data = NPS_DF))
  folds = createFolds(PS_DF$y,10)
  hyper_params = CV_func(PS_DF,w_hat,folds)
  
  ################ MODEL CHUNK #################
  PS_y = PS_DF$y
  NPS_y = NPS_DF$y
  ps_mat  = model.matrix(y ~ .-1, data = PS_DF)
  nps_mat =  model.matrix(y ~ .-1, data = NPS_DF)
  
  #calculate optimal penalty parameters instead of using the cross-validation for beta_hat_optimal
  lambda_star= ((3/ps_n)*popres23) / (5.25*(1-r^2))
  eta_star = lambda_star*r*(sqrt(5.25)/sqrt(w_hat %*% w_hat))
  
  #Fit all estimators on the PS sample (and nps when needed) using the found or calculated penalty parameters 
  ABTLE = solve( t(ps_mat) %*% ps_mat + (hyper_params[[1]]*ps_n * diag(ncol(ps_mat))) ) %*% ( t(ps_mat) %*% PS_y + (hyper_params[[2]]*ps_n * w_hat) )
  RRE = solve( t(ps_mat) %*% ps_mat + (hyper_params[[3]]*ps_n * diag(ncol(ps_mat))) ) %*% ( t(ps_mat) %*% PS_y )
  DBTLE = solve( t(ps_mat) %*% ps_mat + (hyper_params[[4]]*ps_n * diag(ncol(ps_mat))) ) %*% ( t(ps_mat) %*% PS_y + (hyper_params[[4]]*ps_n * w_hat) )
  MLE_ps = solve(t(ps_mat) %*% ps_mat) %*% ( t(ps_mat) %*% PS_y)
  w = solve(t(nps_mat) %*% nps_mat) %*% ( t(nps_mat) %*% NPS_y)
  Bayes_cp = gen.res(PS_y,ps_mat,NPS_y,nps_mat)
  ABTLE_theoretical = solve( t(ps_mat) %*% ps_mat + (lambda_star*ps_n * diag(ncol(ps_mat))) ) %*% ( t(ps_mat) %*% PS_y + (eta_star*ps_n * w_hat) )
  
  #also add the calculated penalty parameters to the hyper_params vector
  hyper_params$lambda_star =lambda_star
  hyper_params$eta_star =as.numeric(eta_star)
  ################ MODEL CHUNK END  #################
  ################ ANALYSIS CHUNK END #################
  return(list("ABTLE" = ABTLE, "RRE" = RRE,"DBTLE"=DBTLE, "MLE_ps" = MLE_ps, "w" = w, "Bayes_cp" = Bayes_cp,"ABTLE_theoretical"=ABTLE_theoretical,"hyper_params"=hyper_params))
}

## 1.7 extract analysis function (litterally just does sapply and makes a new list )
extract_beta_values = function(df) {
  ABTLE = sapply(df, function(x) x$ABTLE)
  RRE = sapply(df, function(x) x$RRE)
  MLE_ps = sapply(df, function(x) x$MLE_ps)
  w = sapply(df, function(x) x$w)
  Bayes_cp = sapply(df, function(x) x$Bayes_cp)
  DBTLE = sapply(df, function(x) x$DBTLE)
  ABTLE_theoretical = sapply(df, function(x) x$ABTLE_theoretical)
  hyper_params = sapply(df, function(x) x$hyper_params)
  return(list(ABTLE = ABTLE,RRE = RRE, DBTLE=DBTLE, MLE_ps = MLE_ps, w = w,Bayes_cp = Bayes_cp,ABTLE_theoretical=ABTLE_theoretical,hyper_params =hyper_params ))
}

## 1.8 calculates the ARMSE and MoAB for each estimator FOR ALL ITERATIONS
eval_func = function(results,betas)  {
  iter_rmse = sapply(results, function(method) {
    # Calculate the RMSE for each row using apply

    apply(method, 2, function(column) sqrt(mean((column - betas)^2)))
    
  })
  iter_armse=colMeans(iter_rmse)
  
  iter_moab = sapply(results, function(method) {
    method_mat=as.matrix(method,nrow=length(betas))
    result = sweep(method_mat, 1, betas, "-")
    result2=mean(abs(rowMeans(result)))
    return(result2)
    
  })
  
  df=data.frame(estimation_method=names(iter_armse),armse=iter_armse,moab=iter_moab,row.names = NULL)
  return(df)
}

# 1.9 Calculates the MSE for each estimator FOR EACH ITERATION
MSE_func = function(results,betas){
  iter_mse = sapply(results, function(method) {
    # Calculate the RMSE for each row using apply
    apply(method, 2, function(column) sqrt(mean((column - betas)^2)))
    
  })
  return(iter_mse)
}

## 2. Body function this functions connects all above functions, running them in order and connecting their outputs It a row out of the scenario df, fixed parameters which never change, and B which determines number of times the analysis is replicated 
simfunction = function(par_comb,par_fixed, B) {
  
  # 1. Import scenario parameters
  p_cov = par_comb$mc
  bias_factor=par_comb$bias_factor
  r=par_comb$r
  y_var=par_comb$y_var
  
  #2. Import fixed scenario parameters
  N = par_fixed$N
  nps_n = par_fixed$nps_n
  p = par_fixed$p
  main_coef=unlist(par_fixed$betas)
  
  # 2. Simulate the population into one big df
  Pop_df = population_func(N, p, p_cov)
  
  # 3. generate phis and add them to a coefs vector with the original coefficients 
  phis=adj_bias_func(main_coef,r,bias_factor,p)
  coefs=c(main_coef,phis)

  # 3. add the population outcome to the population df
  Pop_df = y_func(Pop_df, coefs,y_var)
  
  # 4. Transform selectivity variable z into p
  Pop_df = Pop_df %>%
    mutate(p = Selectivity_func(z, nps_n,2,1)) 
  
  # 4.1 fit the model on the population to get the true residual variance which is needed for the theoretical version of ABTLE 
  pop_model=lm(y~ -1 +X1 + X2+ X3,data=Pop_df)
  popresvar=var(resid(pop_model))
  
  # 5. conduct the analysis B times for each of the four sample sizes 
  ps_n_values = list(25, 50, 75, 100)
  results_list = lapply(ps_n_values, function(ps_n) {
    as.data.frame(replicate(B, analysis_func(Pop_df, "p", ps_n, nps_n, main_coef,r,popresvar)))
  })
  
  #prepare data frames for each of the three types of outputs
  combined_df = data.frame()
  MSE_df = data.frame()
  hyper_comb_df = data.frame()
  #exact and put each output in its correct df (basically just sorting differently sized results into DFs)
  for (i in 1:length(results_list)){
    results_permethod = extract_beta_values(results_list[[i]])
    hyper_out = results_permethod[8]

    #Evaluate each method
    results_permethod=results_permethod[1:7]
    output=eval_func(results_permethod,main_coef)
    combined_df = bind_rows(combined_df, output)
    
    
    hyper_params_ouput=data.frame("lambda_ab" = unlist(hyper_out$hyper_params[1,]),"eta_ab" = unlist(hyper_out$hyper_params[2,]),"lambda_RRE"= unlist(hyper_out$hyper_params[3,]),"lambda_db"= unlist(hyper_out$hyper_params[4,]),"lambda_star"=unlist(hyper_out$hyper_params[5,]),"eta_star"=unlist(hyper_out$hyper_params[6,]) )
    hyper_comb_df = bind_rows(hyper_comb_df, hyper_params_ouput)
    
    MSE_output=as.data.frame(MSE_func(results_permethod,main_coef))
    MSE_df=bind_rows(MSE_df,MSE_output)
    
  }
  combined_df$ps_n=rep(c(25,50,75,100),each=7)
  combined_df$coef=list(coefs)
  #return three df 1, Results relevant for the analysis 2.underlying penalty paramters 3. MSE per iteration
  return(list("combined_df"=combined_df,"hyper_comb_df"=hyper_comb_df, "MSE_df"=MSE_df))
}