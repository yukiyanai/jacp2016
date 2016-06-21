## multilevelmodel.R
########################################################
## To replicate the results reported in
## Yanai's paper for JACP 2016
## Section 4.2: Heterogeneous Treatment Effects
##
## 2016-06-21 Yuki Yanai
########################################################

## Load (or install and load) some packages
if (!require('devtools')) {
    install.packages('devtools', dependencies = TRUE)
    library('devtools')
}
if (!require('dplyr')) {
    install.pakcages('dplyr', dependencies = TRUE)
    library('dplyr')
}
if (!require('rethinking')) {
    install_github('rmcelreath/rethinking')
    library('rethinking')
}

## Create directory to save figures
dir.create('figs', showWarnings = FALSE)

## Set parameter values for the simulation
a <- -0.5
sd_alpha <- 1
b <- 0.8
sd_beta <- 0.2
n_countries <- 60
respondents <- as.integer(rep(c(10, 50, 100), each = 20))

## Data generation
set.seed(2016-06-25)  # use this seed if you want the exactly same result reported
alpha <- rnorm(n_countries, mean = a, sd = sd_alpha)
beta <- rnorm(n_countries, mean = b, sd = sd_beta)
betadf <- data.frame(beta)
quartz(file = 'figs/sim1-true-beta.pdf', type = 'pdf',
       height = 3, width = 4.5)
par(mgp = c(2.2, 1, 0), mar = c(4, 4, 2, 2), cex = .8)
betaplot <- ggplot(betadf, aes(beta)) + geom_histogram(color = "black", binwidth = .1) +
    xlab(expression(beta[j]))
print(betaplot)
dev.off()

country_data <- vector('list', n_countries)
for (i in 1:n_countries) {
    n_res <- respondents[i]
    treatment <- sample(rep(c(0, 1), each = n_res/2), n_res)
    prob <- logistic(alpha[i] + beta[i] * treatment)
    outcome <- rbinom(n_res, size = 1, prob = prob)
    country_data[[i]] <- list(country = i,
                              n_res = n_res,
                              outcome = outcome,
                              treatment = treatment,
                              true_alpha = alpha[i],
                              true_beta = beta[i],
                              true_prob = prob)
}
df <- do.call('rbind', lapply(country_data, data.frame))
# glimpse(df)
# write.csv(df, file = 'jacp2016-simulation1.csv', row.names = FALSE)

## Complete pooling 
glm_pool <- glm(outcome ~ treatment, data = df, family = binomial('logit'))
summary(glm_pool)  

## No pooling (separating)
glm_sep <- glm(outcome ~ treatment * as.factor(country), data = df,
               family = binomial(link = 'logit'))
summary(glm_sep)
a_est_glm_sep <- coef(glm_sep)[1] + c(0, coef(glm_sep)[3:61])
b_est_glm_sep <- coef(glm_sep)[2] + c(0, coef(glm_sep)[62:120])
# weighted mean of estimates
sum(b_est_glm_sep * respondents) / nrow(df)

## Country-fixed effect
glm_fe <- glm(outcome ~ treatment + as.factor(country), data = df,
              family = binomial(link = 'logit'))
summary(glm_fe)
a_est_glm_fe <- coef(glm_fe)[1] + c(0, coef(glm_fe)[3:61])
b_est_glm_fe <- coef(glm_fe)[2]

## Multilevel Model (Partial pooling)
fit_part <- map2stan(
    alist(
        outcome ~ dbinom(1, theta),
        logit(theta) <- alpha[country] + beta[country] * treatment,
        alpha[country] ~ dnorm(a, sa),
        a ~ dnorm(0, 1),
        sa ~ dcauchy(0, 1),
        beta[country] ~ dnorm(b, sb),
        b ~ dnorm(0, 1),
        sb ~ dcauchy(0, 1)
    ),
    data = df, iter = 4000, warmup = 1000,
    chains = 4, cores = 4,
    control = list(adapt_delta = 0.95)
)
## To re-run with new data
#fit_part <- map2stan(fit_part, data = newdata, 
#                    iter = 4000, warmup = 1000,
#                    chains = 4, cores = 4, 
#                    control = list(adapt_delta = 0.95))

## Compare the estimated coefficients between no-pooling and partial pooling
cbind(coef(fit_part), c(coef(fit_sep)[1:60], rep(NA, 2), coef(fit_sep)[61:120], rep(NA, 2)))

## Calculate estiamted probability of success
# Complete pooling
df$p_pool = logistic(coef(glm_pool)[1] + coef(glm_pool)[2] * treatment)
# No-pooling (separating)
df$p_sep <- logistic(a_est_glm_sep[df$country] + b_est_glm_sep[df$country] * df$treatment)
# Fixed effect
df$p_fe <- logistic(a_est_glm_fe[df$country] + b_est_glm_fe * df$treatment)
# Partial poolin (MLM)
est_alpha_part <- coef(fit_part)[1:60]
est_beta_part <- coef(fit_part)[63:122]
## Mean treatment effect
sum(est_beta_part * respondents) / nrow(df)
df$p_part <- logistic(est_alpha_part[df$country] + est_beta_part[df$country] * df$treatment)

## Eaculate absolute error of the estimated probabilities
df$error_pool <- abs(df$p_pool - df$true_prob)
df$error_sep <- abs(df$p_sep - df$true_prob)
df$error_fe <- abs(df$p_fe - df$true_prob)
df$error_part <- abs(df$p_part - df$true_prob)
## Mean error for each country
abs_error <- df %>%
    group_by(country) %>%
    summarize(pooling = mean(error_pool),
              separating = mean(error_sep),
              fixed = mean(error_fe),
              partial = mean(error_part))
## Mean error by the n of observations within country
error_by_size <- df %>% 
    group_by(as.factor(n_res)) %>%
    summarize(pooling = mean(error_pool),
              separating = mean(error_sep),
              fixed = mean(error_fe),
              partial = mean(error_part))
## Plot of absolute error
quartz(file = 'figs/sim1-absolute-error.pdf', type = 'pdf',
       height = 3, width = 4.5)
par(mgp = c(2.2, 1, 0), mar = c(4, 4, 2, 2), cex = .8)
plot(pooling ~ country, data = abs_error,
     xlab = 'country', ylab = 'absolute error', pch = 16,
     ylim = c(0, max(abs_error$pooling)),
     xaxt = 'n')
abline(h = 0, col = 'gray')
axis(1, at = c(1, 20, 40, 60))
points(abs_error$country, abs_error$separating, col = 'royalblue', pch = 16)
points(abs_error$country, abs_error$fixed, col = 'tomato', pch = 16)
points(abs_error$country, abs_error$partial)
text(10, 0.45, 'small n (10)')
text(30, 0.45, 'medium n (50)')
text(50, 0.45, 'large n (100)')
abline(v = c(20.5, 40.5), lwd = .5)
lines(x = c(1, 20), y = rep(error_by_size$pooling[1], 2), col = 'black')
lines(x = c(1, 20), y = rep(error_by_size$separating[1], 2), col = 'royalblue')
lines(x = c(1, 20), y = rep(error_by_size$fixed[1], 2), col = 'tomato')
lines(x = c(1, 20), y = rep(error_by_size$partial[1], 2), lty = 2)
lines(x = c(21, 40), y = rep(error_by_size$pooling[2], 2), col = 'black')
lines(x = c(21, 40), y = rep(error_by_size$separating[2], 2), col = 'royalblue')
lines(x = c(21, 40), y = rep(error_by_size$fixed[2], 2), col = 'tomato')
lines(x = c(21, 40), y = rep(error_by_size$partial[2], 2), lty = 2)
lines(x = c(41, 60), y = rep(error_by_size$pooling[3], 2), col = 'black')
lines(x = c(41, 60), y = rep(error_by_size$separating[3], 2), col = 'royalblue')
lines(x = c(41, 60), y = rep(error_by_size$fixed[3], 2), col = 'tomato')
lines(x = c(41, 60), y = rep(error_by_size$partial[3], 2), lty = 2)
dev.off()

## Data by country
data_by_country <- df %>%
    group_by(country) %>%
    summarize(true_alpha = mean(true_alpha),
              true_beta = mean(true_beta),
              true_prob = mean(true_prob))
data_by_country$alpha_est_pool <- rep(coef(glm_pool)[1], 60)
data_by_country$alpha_est_sep <- a_est_glm_sep
data_by_country$alpha_est_fe <- a_est_glm_fe
data_by_country$alpha_est_part <- coef(fit_part)[1:60]
data_by_country$beta_est_pool <- rep(coef(glm_pool)[2], 60)
data_by_country$beta_est_sep <- b_est_glm_sep
data_by_country$beta_est_fe <- rep(b_est_glm_fe, 60)
data_by_country$beta_est_part <- coef(fit_part)[63:122]


## Treatment effects by country
quartz(file = 'figs/sim1-te-1.pdf', type = 'pdf',
       height = 3, width = 4.5)
par(mgp = c(2.2, 1, 0), mar = c(4, 4, 2, 2), cex = .8)
plot(true_beta ~ country, data = data_by_country,
     xlab = 'country', ylab = 'estimated treatment effect', pch = 16, col = "black",
     ylim = c(-20, 20),
     xaxt = 'n')
axis(1, at = c(1, 20, 40, 60))
abline(h = coef(glm_pool)[2], lty = 2)
abline(h = b_est_glm_fe, col = 'tomato')
points(beta_est_sep ~ country, data = data_by_country, col = 'royalblue', pch = 16)
points(beta_est_part ~ country, data = data_by_country)
text(10, -18, 'small n (10)')
text(30, -18, 'medium n (50)')
text(50, -18, 'large n (100)')
abline(v = c(20.5, 40.5), lwd = .5)
dev.off()
# narrower range
quartz(file = 'figs/sim1-te-2.pdf', type = 'pdf',
       height = 3, width = 4.5)
par(mgp = c(2.2, 1, 0), mar = c(4, 4, 2, 2), cex = .8)
plot(true_beta ~ country, data = data_by_country,
     xlab = 'country', ylab = 'estimated treatment effect', pch = 16, col = "black",
     ylim = c(-.5, 2.5),
     xaxt = 'n')
axis(1, at = c(1, 20, 40, 60))
abline(h = coef(glm_pool)[2], lty = 2)
abline(h = b_est_glm_fe, col = 'tomato')
points(beta_est_sep ~ country, data = data_by_country, col = 'royalblue', pch = 16)
points(beta_est_part ~ country, data = data_by_country)
text(10, 2.4, 'small n (10)')
text(30, 2.4, 'medium n (50)')
text(50, 2.4, 'large n (100)')
abline(v = c(20.5, 40.5), lwd = .5)
dev.off()


## Intercepts by country
quartz(file = 'figs/sim1-intercept-1.pdf', type = 'pdf',
       height = 3, width = 4.5)
par(mgp = c(2.2, 1, 0), mar = c(4, 4, 2, 2), cex = .8)
plot(true_alpha ~ country, data = data_by_country,
     xlab = 'country', ylab = 'intercept', pch = 16, col = "black",
     ylim = c(-20, 20),
     xaxt = 'n')
axis(1, at = c(1, 20, 40, 60))
abline(h = coef(glm_pool)[1])
points(alpha_est_sep ~ country, data = data_by_country, col = 'royalblue', pch = 16)
points(alpha_est_fe ~ country, data = data_by_country, col = 'tomato', pch = 16)
points(alpha_est_part ~ country, data = data_by_country)
text(10, 18, 'small n (10)')
text(30, 18, 'medium n (50)')
text(50, 18, 'large n (100)')
abline(v = c(20.5, 40.5), lwd = .5)
dev.off()
# narrower range
quartz(file = 'figs/sim1-intercept-2.pdf', type = 'pdf',
       height = 3, width = 4.5)
par(mgp = c(2.2, 1, 0), mar = c(4, 4, 2, 2), cex = .8)
plot(true_alpha ~ country, data = data_by_country,
     xlab = 'country', ylab = 'intercept', pch = 16, col = "black",
     ylim = c(-2, 2),
     xaxt = 'n')
axis(1, at = c(1, 20, 40, 60))
abline(h = coef(glm_pool)[1])
points(alpha_est_sep ~ country, data = data_by_country, col = 'royalblue', pch = 16)
points(alpha_est_fe ~ country, data = data_by_country, col = 'tomato', pch = 16)
points(alpha_est_part ~ country, data = data_by_country)
text(10, 1.8, 'small n (10)')
text(30, 1.8, 'medium n (50)')
text(50, 1.8, 'large n (100)')
abline(v = c(20.5, 40.5), lwd = .5)
dev.off()

