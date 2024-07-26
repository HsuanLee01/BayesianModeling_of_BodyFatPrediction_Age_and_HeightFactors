# Bayes Assignment
setwd("C:/Users/Hsuan Lee/Desktop/Methodology and Statistics/Semester 2/Bayesian statistics/Assignment")
options(scipen=999)
data <- read.csv("Body_Fat_Prediction.csv")
data[,1:3] <- scale(data[,1:3], scale = F)

#................
# Gibbs Samplers
#................
Gibbs_Sampler <- function(data, inits, n.chains=2, burn, n.iter,
                          # specify prior
                          mu00 = 0, sigma_square00 = 100000, # intercept prior(normal dist)
                          mu10 = 0, sigma_square10 = 100000, # beta1 prior(normal dist)
                          mu20 = 0, sigma_square20 = 100000, # beta2 prior(normal dist)
                          Alpha0 = 0.001, Beta0 = 0.001 # residual variance prior(gamma dist)
                          ){
  # assign an empty list to each chain
  chain_list <- vector(mode='list', length = n.chains)
  for (chain in 1:n.chains) {
    # assign an empty matrix for sampling in each chain
    keepers <- matrix(0, nrow = n.iter + burn, ncol = 4)
    colnames(keepers)<-c("beta0","beta1","beta2","sigma_square")
    # initial values to different chains
    beta0 <- inits[[chain]][1]
    beta1 <- inits[[chain]][2]
    beta2 <- inits[[chain]][3]
    sigma_square <- inits[[chain]][4]
    # start sampling from the conditional posterior:
    for (iter in 1:(n.iter + burn)) {
      # posterior mean and variance of beta0
      mu01 <- (sum(data[,1] - beta1*data[,2] - beta2*data[,3])/sigma_square + mu00/sigma_square00) /
        (nrow(data)/sigma_square + 1/sigma_square00)
      sigma_square01 <- 1/(nrow(data)/sigma_square + 1/sigma_square00)
      # sample beta0 from the posterior distribution
      beta0 <- rnorm(1, mean = mu01, sd = sqrt(sigma_square01))
      
      # posterior mean and variance of beta1
      mu11 <- (sum(data[,2]*(data[,1] - beta0 - beta2*data[,3])) / sigma_square + (mu10 / sigma_square10)) /
        ((sum((data[,2])^2) / sigma_square) + (1 / sigma_square10))
      sigma_square11 <- 1 / (((sum((data[,2])^2)) / sigma_square) + (1 / sigma_square10))
      # sample beta0 from the posterior distribution
      beta1 <- rnorm(1, mean = mu11, sd = sqrt(sigma_square11))
      
      # posterior mean and variance of beta2
      mu21 <- (sum(data[,3]*(data[,1] - beta0 - beta1*data[,2])) / sigma_square + (mu20 / sigma_square20)) /
        ((sum((data[,3])^2) / sigma_square) + (1 / sigma_square20))
      sigma_square21 <- 1 / (((sum((data[,3])^2)) / sigma_square) + (1 / sigma_square20))
      # sample beta0 from the posterior distribution
      beta2 <- rnorm(1, mean = mu21, sd = sqrt(sigma_square21))
      
      # posterior of sigma square
      Alpha1 <- nrow(data)/2 + Alpha0
      Beta1 <- sum((data[,1] - (beta0 + beta1*data[,2] + beta2*data[,3]))^2) / 2 + Beta0
      # sample sigma square from the posterior distribution
      sigma_square <- MCMCpack::rinvgamma(1, shape = Alpha1, scale = Beta1)
      
      # store the sampled values into keepers
      keepers[iter,] <- c(beta0, beta1, beta2, sqrt(sigma_square))
    }
    # remove the burn in period
    keepers <- keepers[-c(1:burn),]
    # store the sample values of each chain
    chain_list[[chain]] <- keepers
  }
  # return the results
  return(chain_list)
}
#...............................................................................

#.............................
# Metropolis-Hastings Sampler
#.............................
MH_Sampler <- function(data, inits, n.chains=2, burn, n.iter,
                          # specify prior
                          mu00 = 0, sigma_square00 = 100000, # intercept prior(normal dist)
                          mu10 = 0, sigma_square10 = 100000, nu10 = 0.001, # beta1 prior(Non-standardised t-dist)
                          mu20 = 0, sigma_square20 = 100000, nu20 = 0.001, # beta2 prior(Non-standardised t-dist)
                          Alpha0 = 0.001, Beta0 = 0.001, # residual variance prior(gamma dist)
                          tau_square_b1 = .5, tau_square_b2 = .5 # tuning parameter for the proposal density of beta1 and beta2
){
  # assign an empty list to each chain
  chain_list <- vector(mode='list', length = n.chains)
  for (chain in 1:n.chains) {
    # assign an empty matrix for sampling in each chain
    keepers <- matrix(0, nrow = n.iter + burn, ncol = 4)
    colnames(keepers)<-c("beta0","beta1","beta2","sigma_square")
    # initial values to different chains
    beta0 <- inits[[chain]][1]
    beta1 <- inits[[chain]][2]
    beta2 <- inits[[chain]][3]
    sigma_square <- inits[[chain]][4]
    # start sampling from the conditional posterior:
    for (iter in 1:(n.iter + burn)) {
      # posterior mean and variance of beta0
      mu01 <- (sum(data[,1] - beta1*data[,2] - beta2*data[,3])/sigma_square + mu00/sigma_square00) /
        (nrow(data)/sigma_square + 1/sigma_square00)
      sigma_square01 <- 1/(nrow(data)/sigma_square + 1/sigma_square00)
      # sample beta0 from the posterior distribution
      beta0 <- rnorm(1, mean = mu01, sd = sqrt(sigma_square01))
      
      # compute the log of proportional conditional posterior of beta1
      post_b1 <- function(b1){
        (2*b1*sum(data[,2]*(data[,1] - beta0 - beta2*data[,3])) - b1^2*sum(data[,2]^2)) / (2*sigma_square) -
          (log(1 + ((b1-mu10)^2 / (nu10*sigma_square10))) * ((nu10+1)/2))
      }
      # sample from the proposal distribution
      beta1_cand <- rnorm(1, mean = beta1, sd = tau_square_b1)
      # sample from a uniform distribution
      log_u_b1 <- log(runif(1, min = 0, max = 1))
      # compute the acceptance ratio
      log_r_b1 <- (post_b1(beta1_cand) - post_b1(beta1))
      # decide whether to accept the candidate value
      if(log_u_b1 <= log_r_b1){beta1 <- beta1_cand}
      
      # compute the log of proportional conditional posterior of beta2
      post_b2 <- function(b2){
        (2*b2*sum(data[,3]*(data[,1] - beta0 - beta1*data[,2])) - b2^2*sum(data[,3]^2)) / (2*sigma_square) -
          (log(1 + ((b2-mu20)^2 / (nu20*sigma_square20))) * ((nu20+1)/2))
      }
      # sample from the proposal distribution
      beta2_cand <- rnorm(1, mean = beta2, sd = tau_square_b2)
      # sample from a uniform distribution
      log_u_b2 <- log(runif(1, min = 0, max = 1))
      # compute the acceptance ratio
      log_r_b2 <- (post_b2(beta2_cand) - post_b2(beta2))
      # decide whether to accept the candidate value
      if(log_u_b2 <= log_r_b2){beta2 <- beta2_cand}
      
      # posterior of sigma square
      Alpha1 <- nrow(data)/2 + Alpha0
      Beta1 <- sum((data[,1] - (beta0 + beta1*data[,2] + beta2*data[,3]))^2) / 2 + Beta0
      # sample sigma square from the posterior distribution
      sigma_square <- MCMCpack::rinvgamma(1, shape = Alpha1, scale = Beta1)
      
      # store the sampled values into keepers
      keepers[iter,] <- c(beta0, beta1, beta2, sqrt(sigma_square))
    }
    # remove the burn in period
    keepers <- keepers[-c(1:burn),]
    # store the sample values of each chain
    chain_list[[chain]] <- keepers
  }
  # return the results
  return(chain_list)
}
#...............................................................................

#.................................
# Assess convergence of the model
#.................................

# 1. Trace Plot

Trace_Plot <- function(sample){
  par(mfrow = c(1, 4))
  for (parameter in 1:4) {
    # chain 1
    plot(sample[[1]][,parameter], type = "l", main = "Trace Plot", xlab = "Iterations", 
         ylab = paste("", colnames(sample[[1]])[parameter]), col = "red")
    # chain 2
    lines(sample[[2]][,parameter], type = "l", col = "blue")
  }
}

#...............................................................................

# 2. Density Plot

Density_Plot <- function(sample){
  par(mfrow = c(1, 4))
  for (parameter in 1:4) {
    plot(density(sample[[1]][,parameter]), main = "Density Plot", 
         xlab = paste("", colnames(sample[[1]])[parameter]), ylab = "Density", col = "red")
    lines(density(sample[[2]][,parameter]), col = "blue")
  }
}
#...............................................................................

# 3. Autocorrelation Plot

Autocor_Plot <- function(sample){
  par(mfrow = c(2, 4))
  for (parameter in 1:4){
    acf(sample[[1]][,parameter], type = "correlation", ylab = "Autocorrelation",
        main = paste("Autocorrelation Plot of", colnames(sample[[1]])[parameter], "(chain 1)"))
    acf(sample[[2]][,parameter], type = "correlation", ylab = "Autocorrelation",
        main = paste("Autocorrelation Plot of", colnames(sample[[2]])[parameter], "(chain 2)"))
  }
}
#...............................................................................

# 4. Gelman-Rubin statistic

Gelman_Plot <- function(sample, n_iter){
  par(mfrow = c(1, 4))
  R_hat <- matrix(0, nrow = n_iter-1, ncol = 4) # store Gelman-Rubin Convergence Diagnostic
  for (parameter in 1:4) {
    thetaj_chain1 <- 0 # store chain 1 mean
    thetaj_chain2 <- 0 # store chain 2 mean
    theta <- 0 # store grand mean
    B <- 0 # store between chain variance
    W <- 0 # store within chain variance
    V <- 0 # store total chain variance
    for (iter in 2:n_iter) {
      # chain mean
      thetaj_chain1 <- mean(sample[[1]][1:iter, parameter])
      thetaj_chain2 <- mean(sample[[2]][1:iter, parameter])
      # grand mean across chain
      theta <- (sum(sample[[1]][1:iter, parameter] + sample[[2]][1:iter, parameter])) / (2*iter)
      # between chain variance
      B <- iter/(length(sample)-1) * ((thetaj_chain1-theta)^2 + (thetaj_chain2-theta)^2)
      # within chain variance
      W <- 1/((iter-1)*length(sample)) * sum((sample[[1]][1:iter, parameter]-thetaj_chain1)^2 + 
        (sample[[2]][1:iter, parameter]-thetaj_chain2)^2)
      # total chain varicance
      V <- (iter-1)/iter * W + (1/iter)*B
      # Gelman-Rubin Convergence Diagnostic
      R_hat[iter-1, parameter] <- sqrt(V/W)
    }
    plot(R_hat[,parameter], type = "l", lwd = 2,
         main = paste("Gelman-Rubin Statistics of", colnames(sample[[2]])[parameter]),
         xlab = "iteration", ylab = expression(hat(R)), ylim = c(1, 1.05))
    abline(h = 1, col = "red", lwd = 2)
  }
}
#...............................................................................

#..................
# Summary function
#..................
Bayes_Summary <- function(sample){
  # create a empty matrix to store the mean, sd and MC error
  summ <- matrix(0, nrow = 4, ncol = 3)
  colnames(summ) <- c("mean", "sd", "MC error")
  rownames(summ) <- c("beta0", "beta1", "beta2", "sigma_square")
  # combine two chains
  tot <- rbind(sample[[1]], sample[[2]])
  # calculate the mean
  summ[, 1] <- apply(tot, 2, mean)
  # calculate sd
  summ[, 2] <- apply(tot, 2, sd)
  # calculate MC error should < .05
  summ[, 3] <- apply(tot, 2, function(x) sd(x)/(2*nrow(sample[[1]])))
  # calculate the credible intervals
  CI <- apply(tot, 2, quantile, prob = c(0.025, 0.25, 0.5, 0.75, 0.975))
  list <- list(Summary = summ, Credible_Interval = CI)
  return(list)
}
#...............................................................................

#.....
# DIC
#.....
DIC <- function(sample, data){
  # combine two chains
  tot <- rbind(sample[[1]], sample[[2]])
  # compute Dhat
  Dhat <- -2 * sum(dnorm(data[,1], 
                         mean = mean(tot[,1]) + mean(tot[,2])*data[,2] + mean(tot[,3])*data[,3],
                         sd = mean(tot[,4]), log = T))
  # compute Dbar
  loglike <- matrix(0, nrow = nrow(tot), ncol = 1)
  for (q in 1:nrow(tot)) {
    loglike[q,] <- sum(dnorm(data[,1],
                         mean = tot[q,1] + tot[q,2]*data[,2] + tot[q,3]*data[,3],
                         sd = tot[q,4], log = T))
  }
  Dbar <- -2 * apply(loglike, 2, mean)
  # compute DIC
  DIC <- 2*Dbar - Dhat
  # estimate of the effective number of parameters
  PD <- Dbar - Dhat
  list <- list(DIC = DIC, PD = PD)
  return(list)
}
#...............................................................................

#..............
# Bayes Factor
#..............
Bayes_Factor <- function(sample, data){
  # extract the estimates of beta1 and beta2
  beta1_est <- Bayes_Summary(sample)[[1]][2,1]
  beta2_est <- Bayes_Summary(sample)[[1]][2,1]
  # standardized the estimates of beta1 and beta2
  beta1_hat <- beta1_est*sd(data[, 2])/sd(data[, 1])
  beta2_hat <- beta2_est*sd(data[, 3])/sd(data[, 1])
  beta_hat <- c(beta1_hat, beta2_hat)
  # compute the covariance matrix
  tot <- rbind(sample[[1]], sample[[2]]) # combine two chains
  sigma_beta12 <- cov(tot[,c("beta1", "beta2")])
  # sample from the posterior
  library(mvtnorm)
  g_sample <- rmvnorm(100000, mean = beta_hat, sigma = sigma_beta12)
  
  # fractional prior distribution
  J <- 2 # number of independent constraints in the hypotheses under consideration
  N <- nrow(data)
  b <- J / N # fraction of the density of the data used to specify a prior distribution
  # sample from the fractional prior distribution
  h_sample <- rmvnorm(100000, mean = c(0, 0), sigma = sigma_beta12 / b)
  
  # evaluate the informative hypotheses
  # H1: |b1 - b2| < .1
  BF_1u <- mean(abs(g_sample[, 1] - g_sample[, 2]) < 0.1) / mean(abs(h_sample[, 1] - h_sample[, 2]) < 0.1)
  # H2: b1 > 0, b2 > 0
  BF_2u <- mean(g_sample[, 1] > 0 & g_sample[, 2] > 0) / mean(h_sample[, 1] > 0 & h_sample[, 2] > 0)
  # Hu: b1, b2
  BF_21 <- BF_2u / BF_1u
  return(list(BF_1u = BF_1u, BF_2u = BF_2u, BF_21 = BF_21))
}
#...............................................................................

##############
## Analysis ##
##############

# specify the initial values
#............................
init1 <- c(beta0 = 0.9, beta1 = 0.2, beta2 = 0.5, sigma_square = 2.0)
init2 <- c(beta0 = 1.8, beta1 = 2.6, beta2 = 1.5, sigma_square = 0.6)
init <- list(init1, init2)

# Gibbs Sampler
#..............
# Implement the Gibbs Sampler to obtain the samples
set.seed(9252568)
samples_gibbs <- Gibbs_Sampler(data = data, inits = init, 
                               n.chains=2, burn = 500, n.iter = 10000)

# inspect the convergence of the samples sampled by the Gibbs Sampler
#Trace_Plot(samples_gibbs)
#Density_Plot(samples_gibbs)
#Autocor_Plot(samples_gibbs)
#Gelman_Plot(samples_gibbs, n_iter = 1500)

# obtain the result of the samples sampled by the Gibbs Sampler
#Bayes_Summary(samples_gibbs)

# However, the samples sampled by the Gibbs Sampler would not be used in this study.
# Since the "non-standardized t distribution" used in MH sampler allows for more extreme values,
# which is a more proper distribution for the study.
#...............................................................................

# MH Sampler Full Model
#......................
# Implement the MH Sampler to obtain the samples
set.seed(9252568)
samples_MH <- MH_Sampler(data = data,inits = init,
                         n.chains = 2, burn = 500, n.iter = 10000)

# inspect the convergence of the samples sampled by the MH Sampler
Trace_Plot(samples_MH)
Density_Plot(samples_MH)
Autocor_Plot(samples_MH)
Gelman_Plot(samples_MH, n_iter = 10000)

# obtain the result of the samples sampled by the MH Sampler
Bayes_Summary(samples_MH)

# obtain DIC and Bayes Factor of the model 
DIC(samples_MH, data = data) # DIC = 1769.945, PD = 3.913394
set.seed(9252568)
Bayes_Factor(samples_MH, data)

# the samples sampled by the MH Sampler would be used in this study.
#...............................................................................

# Intercept only model
#......................................
# create the data without variable "Age" and "Height"
data_emp <- data
data_emp[,c(2,3)] <- 0

# Implement the MH Sampler to obtain the samples
set.seed(9252568)
samples_emp <- MH_Sampler(data = data_emp,inits = init,
                          n.chains = 2, burn = 500, n.iter = 10000)

# obtain DIC of the model 
DIC(samples_emp, data = data_emp) # DIC = 1788.898, PD = 1.996333

#...............................................................................

# The model without variable "Height"
#......................................
# create the data without variable "Height"
data_age <- data
data_age[,3] <- 0

# Implement the MH Sampler to obtain the samples
set.seed(9252568)
samples_age <- MH_Sampler(data = data_age,inits = init,
                          n.chains = 2, burn = 500, n.iter = 10000)

# obtain DIC of the model 
DIC(samples_age, data = data_age) # DIC = 1768.376, PD = 2.915565

# obtain the result of the samples sampled by the MH Sampler
Bayes_Summary(samples_age)
#...............................................................................

# Frequentist Linear Regression Model
#.....................................
LRM <- lm(BodyFat~ Age + Height, data = data)
summary(LRM)
