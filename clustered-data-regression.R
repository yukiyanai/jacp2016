## clustered-data-regression.R
#####################################################
## Replicate the results reported in
## Yanai's paper for JACP 2016
## Section 4.1 Regression with Clustered Data
##
## 2016-06-21 Yuki Yanai
#####################################################

## Load some packages. Install before loading if necessary.
library('arm')
library('mvtnorm')
library('lme4')
library('multiwayvcov')
library('ggplot2')
library('dplyr')
library('haven')

## Create directory to save figures
dir.create('figs', showWarnings = FALSE)

## Setting Up Simulations

## Prepare some functions
gen_cluster <- function(param = c(.1, .5), n = 1000, n_cluster = 50, rho = .5) {
    # Function to generate clustered data
    # Required package: mvtnorm
    
    # individual level
    Sigma_i <- matrix(c(1, 0, 0, 1 - rho), ncol = 2)
    values_i <- rmvnorm(n = n, sigma = Sigma_i)
    
    # cluster level
    cluster_name <- rep(1:n_cluster, each = n / n_cluster)
    Sigma_cl <- matrix(c(1, 0, 0, rho), ncol = 2)
    values_cl <- rmvnorm(n = n_cluster, sigma = Sigma_cl)
    
    # predictor var consists of individual- and cluster-level components
    x <- values_i[ , 1] + rep(values_cl[ , 1], each = n / n_cluster)
    
    # error consists of individual- and cluster-level components
    error <- values_i[ , 2] + rep(values_cl[ , 2], each = n / n_cluster)
    
    # data generating process
    y <- param[1] + param[2]*x + error
    
    df <- data.frame(x, y, cluster = cluster_name)
    return(df)
}

## Generate a dataset with clusters and fit OLS
# Calculate cluster-robust SE when cluster_robust = TRUE
cluster_sim <- function(param = c(.1, .5), n = 1000, n_cluster = 50,
                            rho = .5, cluster_robust = FALSE) {
    # Required packages: mvtnorm, multiwayvcov
    df <- gen_cluster(param = param, n = n , n_cluster = n_cluster, rho = rho)
    fit <- lm(y ~ x, data = df)
    b1 <- coef(fit)[2]
    if (!cluster_robust) {
        Sigma <- vcov(fit)
        se <- sqrt(diag(Sigma)[2])
        b1_ci95 <- confint(fit)[2, ]
    } else { # cluster-robust SE
        Sigma <- cluster.vcov(fit, ~ cluster)
        se <- sqrt(diag(Sigma)[2])
        t_critical <- qt(.025, df = n - 2, lower.tail = FALSE)
        lower <- b1 - t_critical*se
        upper <- b1 + t_critical*se
        b1_ci95 <- c(lower, upper)
    }
    return(c(b1, se, b1_ci95))
}

# Function to iterate the simulation. A data frame is returned.
run_cluster_sim <- function(n_sims = 1000, param = c(.1, .5), n = 1000,
                            n_cluster = 50, rho = .5, cluster_robust = FALSE) {
    # Required packages: mvtnorm, multiwayvcov, dplyr
    df <- replicate(n_sims, cluster_sim(param = param, n = n, rho = rho,
                                        n_cluster = n_cluster,
                                        cluster_robust = cluster_robust))
    df <- as.data.frame(t(df))
    names(df) <- c('b1', 'se_b1', 'ci95_lower', 'ci95_upper')
    df <- df %>% 
        mutate(id = 1:n(),
               param_caught = ci95_lower <= param[2] & ci95_upper >= param[2])
    return(df)
}

## Genrate data without Clusters
sim_params <- c(.4, 0)   # beta1 = 0: no effect of x on y
sim_nocluster <- run_cluster_sim(n_sims = 10000, param = sim_params, rho = 0)
hist_nocluster <- ggplot(sim_nocluster, aes(b1)) +
    geom_histogram(color = 'black') +
    geom_vline(xintercept = sim_params[2], color = 'red')
print(hist_nocluster)
## Confidence intervals
ci95_nocluster <- ggplot(sample_n(sim_nocluster, 100),
                         aes(x = reorder(id, b1), y = b1, 
                             ymin = ci95_lower, ymax = ci95_upper,
                             color = param_caught)) +
    geom_hline(yintercept = sim_params[2], linetype = 'dashed') +
    geom_pointrange() +
    labs(x = 'sim ID', y = 'b1', title = 'Randomly Chosen 100 95% CIs') +
    scale_color_discrete(name = 'True param value', labels = c('missed', 'hit')) +
    coord_flip()
print(ci95_nocluster)
## Rate of type I error
sim_nocluster %>% summarize(type1_error = 1 - sum(param_caught)/n())

## Generate data with clusters
set.seed(2016-06-25)  # use this seed to get the exactly same results
sim_params <- c(.4, 0)   # beta1 = 0: no effect of x on y
sim_cluster_ols <- run_cluster_sim(n_sims = 10000, param = sim_params)
## Histogram of parameter estiamtes
hist_cluster_ols <- hist_nocluster %+% sim_cluster_ols +
    xlab(expression(paste('estimated ', beta)))
quartz(file = 'figs/cluster-hist-unbiased.pdf', type = 'pdf',
       height = 3, width = 3)
print(hist_cluster_ols)
dev.off()
## Histogram of standard error
hist_se <- ggplot(sim_cluster_ols, aes(se_b1)) + 
    geom_histogram(color = "black") +
    xlab(expression(paste('se of ', beta)))
quartz(file = 'figs/cluster-hist-se.pdf', type = 'pdf',
       height = 3, width = 3)
print(hist_se)
dev.off()
## 95% CI
ci95_cluster_ols <- ci95_nocluster %+% sample_n(sim_cluster_ols, 100)
print(ci95_cluster_ols)
## Rate of type I error
sim_cluster_ols %>% summarize(type1_error = 1 - sum(param_caught)/n())


## Cluster-robust SE
sim_params <- c(.4, 0)   # beta1 = 0: no effect of x on y
sim_cluster_robust <- run_cluster_sim(n_sims = 10000, param = sim_params,
                                      cluster_robust = TRUE)
## Histogram of coefficient estimate
hist_cluster_robust <- hist_nocluster %+% sim_cluster_ols
print(hist_cluster_robust)
## Histogram of cluster-robust SE
hist_crse <- ggplot(sim_cluster_robust, aes(se_b1)) + 
    geom_histogram(color = "black") +
    xlab(expression(paste('se of ', beta)))
quartz(file = 'figs/cluster-hist-crse.pdf', type = 'pdf',
       height = 3, width = 3)
print(hist_crse)
dev.off()
## 95% CI with cluster-robust SE
ci95_cluster_robust <- ci95_nocluster %+% sample_n(sim_cluster_robust, 100)
print(ci95_cluster_robust)
## type I error
sim_cluster_robust %>% summarize(type1_error = 1 - sum(param_caught)/n())


## Cluster-Robust SE, Fixed Effect, or Multilevel Models

## Function to compare 4 methods
sim_4models <- function(param = c(.1, .5), n = 1000, rho = .5, n_cluster = 50) {
    # Required packages: mvtnorm, multiwaycov, lme4,
    df <- gen_cluster(param = param, n = n , n_cluster = n_cluster, rho = rho)
    
    # Fit three different models
    fit_ols <- lm(y ~ x, data = df)
    fit_fe <- lm(y ~ x + factor(cluster), data = df)
    fit_multi <- lmer(y ~ x + (1|cluster), data = df)
    
    coef <- c(coef(fit_ols)[2], coef(fit_fe)[2], fixef(fit_multi)[2])
    
    # 4 ways to get var-cov matrix
    Sigma_ols <- vcov(fit_ols)
    Sigma_crse <- cluster.vcov(fit_ols, ~ cluster)
    Sigma_fe <- vcov(fit_fe)
    Sigma_multi <- vcov(fit_multi)
    
    # 4 different SEs
    se <- sqrt(c(diag(Sigma_ols)[2], diag(Sigma_crse)[2],
                 diag(Sigma_fe)[2], diag(Sigma_multi)[2]))
    
    return(c(coef[1], se[1],  # OLS
             coef[1], se[2],  # OLS with cluster-robust SE
             coef[2], se[3],  # fixed-effect model
             coef[3], se[4])  # multi-level model
           )
}

## Run the simulation
set.seed(2016-06-24)  # use this seed to get the exactly same result
n_sims <- 10000
sim_params <- c(.4, 0)
comp_4models <- replicate(n_sims, sim_4models(param = sim_params))
comp_df <- as.vector(comp_4models)
comp_df <- matrix(comp_df, ncol = 2, byrow = TRUE)
comp_df <- as.data.frame(comp_df)
names(comp_df) <- c('b1', 'se_b1')
comp_df <- comp_df %>%
    mutate(model = rep(c('OLS', 'CRSE', 'FE', 'MLM'), n_sims),
           id = rep(1:n_sims, each = 4))

## Compare coefficient estimates
density_b1 <- ggplot(comp_df, aes(b1, color = model)) +
    geom_density() +
    xlab(expression(paste('estimate of ', beta)))
quartz(file = 'figs/compare-4models.pdf', type = 'pdf',
       height = 3, width = 4.5)
print(density_b1)
dev.off()

## Compare type 1 error rates
comp_df <- comp_df %>%
    mutate(lower = b1 - 1.96*se_b1,
           upper = b1 + 1.96*se_b1,
           param_caught = lower <= sim_params[2] & upper >= sim_params[2])
comp_df %>% group_by(model) %>%
    summarize(type1_error = 1 - sum(param_caught)/n())
