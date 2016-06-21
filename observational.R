## observational.R
############################################################
## To replicate the results reported in
## Yanai's paper for JACP 2016
## Section 4.3: Multilevel Models in Observational Studies
##
## 2016-06-21 Yuki Yanai
############################################################

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
if (!require('multiwayvcov')) {
    install.packages('multiwayvcov', dependencies = TRUE)
    library('multiwayvcov')
}
if (!require('MatchIt')) {
    install.packages('MatchIt', dependencies = TRUE)
    library('MatchIt')
}
if (!require('optmatch')) {
    install.packages('optmatch', dependencies = TRUE)
    library('optmatch')
}
if (!require('pbkrtest')) {
    install.packages('pbkrtest', dependencies = TRUE)
    library('pbkrtest')
}


## Data generation
n_countries <- 60
n_units <- as.integer(rep(c(10, 50, 100), each = 20))

set.seed(2016-06-24) # use this seed to get the exactly same results
te <- rnorm(n_countries, 2, sd = 1)
gamma1 <- rnorm(n_countries, 0, 4)
gamma2 <- rnorm(n_countries, 0, 4)
gamma3 <- rnorm(n_countries, 0, 4)
gamma4 <- rnorm(n_countries, 0, 4)
z1 <- rnorm(n_countries, 0, 1)
z2 <- rnorm(n_countries, 0, 1)
alpha <- rnorm(n_countries, 2, 2)
beta1 <- 1
beta2 <- -1

country_data <- vector('list', n_countries)
for (i in 1:n_countries) {
    n <- n_units[i]
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    p <- logistic(gamma1[i] * x1 + gamma2[i] *x2 + gamma3[i] * z1[i] + gamma4[i] * z2[i])
    propensity <- 0.5 + (p - 0.5)/2
    treatment <- rbinom(n, size = 1, prob = propensity)
    outcome <- alpha[i] + te[i] * treatment + beta1*x1 + beta2*z1[i] + rnorm(n, 0, 1)
    country_data[[i]] <- list(country = rep(i, n),
                              n_units = rep(n, n),
                              outcome = outcome,
                              treatment = treatment,
                              x1 = x1, x2 = x2, 
                              z1 = rep(z1[i], n), z2 = rep(z2[i], n),
                              true_g1 = rep(gamma1[i], n),
                              true_g2 = rep(gamma2[i], n),
                              true_g3 = rep(gamma3[i], n),
                              true_g4 = rep(gamma4[i], n),
                              true_alpha = rep(alpha[i], n),
                              true_b1 = rep(beta1, n),
                              true_b2 = rep(beta2, n),
                              true_te = rep(te[i], n),
                              true_propensity = propensity
                              )
}
df <- do.call('rbind', lapply(country_data, data.frame))
df$torc <- ifelse(df$treatment == 1, "treatment", "control")

## Proportion of treated by country
df %>% group_by(country) %>%
    summarize(treat_prop = mean(treatment)) %>%
    summary()


## Histogram of true treatment effect
df_country <- df %>% 
    group_by(country) %>%
    summarize(true_te = mean(true_te),
              true_alpha = mean(true_alpha),
              n_units = mean(n()),
              prop_treat = mean(treatment))
te_hist <- ggplot(df_country, aes(true_te)) +
    geom_histogram(binwidth = .25, color = "black") +
    xlab(expression(tau[j]))
quartz(file = 'figs/sim2-true-te.pdf', type = 'pdf',
       height = 3, width = 4.5)
print(te_hist)
dev.off()


## Complete pooling
fit_pool <- lm(outcome ~ treatment + x1 + x2 + z1 + z2, data = df)
Sigma <- cluster.vcov(fit_pool, ~ country)
crse_pool <- sqrt(diag(Sigma))
summary(fit_pool)
crse_pool
te_pool <- rep(coef(fit_pool)[2], n_countries)
a_pool <- rep(coef(fit_pool)[1], n_countries)
mean(te_pool)

## Fixed effect model
fit_fe <- lm(outcome ~ treatment + x1 + x2 + z1 + z2 + as.factor(country), data = df)
Sigma <- cluster.vcov(fit_fe, ~ country)
crse_fe <- sqrt(diag(Sigma))
summary(fit_fe)
crse_fe[2]
te_fe <- rep(coef(fit_fe)[2], n_countries)
mean(te_fe)

## Separaing (no pooling)
fit_sep <- lm(outcome ~ treatment*as.factor(country) + x1 + x2 + z1 + z2, data = df)
summary(fit_sep)
re <- c(0, coef(fit_sep)[66:124])
te_sep <- coef(fit_sep)[2] + re
mean(te_sep, na.rm = TRUE)

## Partial pooling: multilevel
fit_part <- map2stan(
    alist(
        outcome ~ dnorm(theta, sigma),
        theta <- alpha[country] + te[country]*treatment + beta1*x1 + beta2*z1 +
            beta3*x2 + beta4*z2,
        alpha[country] ~ dnorm(a1, s1),
        a1 ~ dnorm(0, 1),
        s1 ~ dcauchy(0, 1),
        te[country] ~ dnorm(a2, s2),
        a2 ~ dnorm(0, 1),
        s2 ~ dcauchy(0, 1),
        beta1 ~ dnorm(0, 1),
        beta2 ~ dnorm(0, 1),
        beta3 ~ dnorm(0, 1),
        beta4 ~ dnorm(0, 1),
        sigma ~ dcauchy(0, 1)
    ),
    data = df, warmup = 1000, iter = 4000,
    chains = 4, cores = 4,
    control = list(adapt_delta = 0.95)
)
# To re-run the same model with new data
#newdata <- list(outcome = df$outcome, treatment = df$treatment, 
#                x1 = df$x1, x2 = df$x2, z1 = df$z1, z2 = df$z2,
#                country = df$country)
#fit_part <- map2stan(fit_part, data = newdata,
#                     warmup = 1000, iter = 4000,
#                     chains = 4, cores = 4,
#                     control = list(adapt_delta = 0.95))
summary(fit_part)
te_part <- coef(fit_part)[63:122]
mean(te_part)

## Compare mean estiamtes
c(mean(te), mean(te_pool), mean(te_fe), mean(te_sep, na.rm = TRUE), mean(te_part))

df_country$te_pool <- te_pool
df_country$te_fe <- te_fe
df_country$te_sep <- te_sep
df_country$te_part <- te_part
df_country <- df_country %>%
    mutate(err_pool = abs(te_pool - true_te),
           err_fe = abs(te_fe - true_te),
           err_sep = abs(te_sep - true_te),
           err_part = abs(te_part - true_te))
error_by_size <- df_country %>% 
    group_by(as.factor(n_units)) %>%
    summarize(pooling = mean(err_pool),
              separating = mean(err_sep, na.rm = TRUE),
              fixed = mean(err_fe),
              partial = mean(err_part))

## Absolute error
quartz(file = 'figs/sim2-absolute-error.pdf', type = 'pdf',
       height = 3, width = 4.5)
par(mgp = c(2.2, 1, 0), mar = c(4, 4, 2, 2), cex = .8)
plot(err_pool ~ country, data = df_country,
     xlab = 'country', ylab = 'absolute error', pch = 16,
     ylim = c(0, 3), xaxt = 'n')
abline(h = 0, col = "gray")
axis(1, at = c(1, 20, 40, 60))
points(err_sep ~ country, data = df_country, col = 'royalblue', pch = 16)
points(err_fe ~ country, data = df_country, col = 'tomato', pch = 16)
points(err_part ~ country, data = df_country)
text(10, 2.8, 'small n (10)')
text(30, 2.8, 'medium n (50)')
text(50, 2.8, 'large n (100)')
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


## Treatment effects by country
quartz(file = 'figs/sim2-te.pdf', type = 'pdf',
       height = 3, width = 4.5)
par(mgp = c(2.2, 1, 0), mar = c(4, 4, 2, 2), cex = .8)
plot(true_te ~ country, data = df_country,
     xlab = 'country', ylab = 'estimated treatment effect', pch = 16, col = "black",
     ylim = c(0, 5),
     xaxt = 'n')
axis(1, at = c(1, 20, 40, 60))
abline(h = coef(fit_pool)[2], lty = 2)
abline(h = coef(fit_fe)[2], col = 'tomato')
points(te_sep ~ country, data = df_country, col = 'royalblue', pch = 16)
points(te_part ~ country, data = df_country)
text(10, 4.7, 'small n (10)')
text(30, 4.7, 'medium n (50)')
text(50, 4.7, 'large n (100)')
abline(v = c(20.5, 40.5), lwd = .5)
dev.off()


## Matching
m_full <- matchit(treatment ~ x1 + x2 + z1 + z2, data = df, method = 'full')
summary(m_full)
summary_m <- summary(m_full, standardize = TRUE)
plot(summary_m)
plot(m_full)
plot(m_full, type = 'jitter')
plot(m_full, type = 'hist')
df_match <- match.data(m_full)

## Matching within each group
mdf <- NULL
for (i in 1:n_countries) {
    df_sub <- df %>% filter(country == i)
    m_out <- matchit(treatment ~ x1 + x2, data = df_sub, method = 'nearest')
    m_sub <- match.data(m_out)
    mdf <- rbind(mdf, m_sub)
}
mdf %>% group_by(country) %>%
    summarize(n_units = mean(n())) %>%
    summary()


## Complete pooling with matched data
mfit_pool <- lm(outcome ~ treatment + x1 + z1, data = mdf)
Sigma <- cluster.vcov(mfit_pool, ~ country)
crse_mpool <- sqrt(diag(Sigma))
summary(mfit_pool)
crse_mpool
te_mpool <- rep(coef(mfit_pool)[2], n_countries)
mean(te_mpool)

## Fixed effect model with matched data
mfit_fe <- lm(outcome ~ treatment + x1 + z1 + as.factor(country), data = mdf)
Sigma <- cluster.vcov(mfit_fe, ~ country)
crse_mfe <- sqrt(diag(Sigma))
summary(mfit_fe)
crse_mfe[2]
te_mfe <- rep(coef(mfit_fe)[2], n_countries)
mean(te_fe)

## Separaing (no pooling) with matched data
mfit_sep <- lm(outcome ~ treatment*as.factor(country) + x1 +  z1, data = mdf)
summary(mfit_sep)
re <- c(0, coef(mfit_sep)[64:122])
te_msep <- coef(mfit_sep)[2] + re
mean(te_msep)

## Partial pooling (multilevel model) with matched data
mfit_part <- map2stan(
    alist(
        outcome ~ dnorm(theta, sigma),
        theta <- alpha[country] + te[country]*treatment + beta1*x1 + beta2*z1,
        alpha[country] ~ dnorm(a1, s1),
        a1 ~ dnorm(0, 1),
        s1 ~ dcauchy(0, 1),
        te[country] ~ dnorm(a2, s2),
        a2 ~ dnorm(0, 1),
        s2 ~ dcauchy(0, 1),
        beta1 ~ dnorm(0, 1),
        beta2 ~ dnorm(0, 1),
        sigma ~ dcauchy(0, 1)
    ),
    data = mdf, warmup = 1000, iter = 4000,
    chains = 4, cores = 4,
    control = list(adapt_delta = 0.95)
)
## To re-run the same model with new data
#newdata <- list(outcome = mdf$outcome, treatment = mdf$treatment, 
#                x1 = mdf$x1, x2 = mdf$x2, z1 = mdf$z1, z2 = mdf$z2,
#                country = mdf$country)
#mfit_part <- map2stan(fit_part, data = newdata,
#                     warmup = 1000, iter = 4000,
#                     chains = 4, cores = 4,
#                     control = list(adapt_delta = 0.95))
summary(mfit_part)
te_mpart <- coef(mfit_part)[63:122]
mean(te_mpart)

## Compare mean estiamtes
c(mean(te), mean(te_mpool), mean(te_mfe), mean(te_msep), mean(te_mpart))

df_country$te_mpool <- te_mpool
df_country$te_mfe <- te_mfe
df_country$te_msep <- te_msep
df_country$te_mpart <- te_mpart
df_country <- df_country %>%
    mutate(err_mpool = abs(te_mpool - true_te),
           err_mfe = abs(te_mfe- true_te),
           err_msep = abs(te_msep - true_te),
           err_mpart = abs(te_mpart - true_te))
error_by_size <- df_country %>% 
    group_by(as.factor(n_units)) %>%
    summarize(mpooling = mean(err_mpool),
              mseparating = mean(err_msep),
              mfixed = mean(err_mfe),
              mpartial = mean(err_mpart))

## Absolute error
quartz(file = 'figs/sim3-absolute-error.pdf', type = 'pdf',
       height = 3, width = 4.5)
par(mgp = c(2.2, 1, 0), mar = c(4, 4, 2, 2), cex = .8)
plot(err_mpool ~ country, data = df_country,
     xlab = 'country', ylab = 'absolute error', pch = 16,
     ylim = c(0, 3), xaxt = 'n')
abline(h = 0, col = "gray")
axis(1, at = c(1, 20, 40, 60))
points(err_msep ~ country, data = df_country, col = 'royalblue', pch = 16)
points(err_mfe ~ country, data = df_country, col = 'tomato', pch = 16)
points(err_mpart ~ country, data = df_country)
text(10, 2.8, 'small n (10)')
text(30, 2.8, 'medium n (50)')
text(50, 2.8, 'large n (100)')
abline(v = c(20.5, 40.5), lwd = .5)
lines(x = c(1, 20), y = rep(error_by_size$mpooling[1], 2), col = 'black')
lines(x = c(1, 20), y = rep(error_by_size$mseparating[1], 2), col = 'royalblue')
lines(x = c(1, 20), y = rep(error_by_size$mfixed[1], 2), col = 'tomato')
lines(x = c(1, 20), y = rep(error_by_size$mpartial[1], 2), lty = 2)
lines(x = c(21, 40), y = rep(error_by_size$mpooling[2], 2), col = 'black')
lines(x = c(21, 40), y = rep(error_by_size$mseparating[2], 2), col = 'royalblue')
lines(x = c(21, 40), y = rep(error_by_size$mfixed[2], 2), col = 'tomato')
lines(x = c(21, 40), y = rep(error_by_size$mpartial[2], 2), lty = 2)
lines(x = c(41, 60), y = rep(error_by_size$mpooling[3], 2), col = 'black')
lines(x = c(41, 60), y = rep(error_by_size$mseparating[3], 2), col = 'royalblue')
lines(x = c(41, 60), y = rep(error_by_size$mfixed[3], 2), col = 'tomato')
lines(x = c(41, 60), y = rep(error_by_size$mpartial[3], 2), lty = 2)
dev.off()


## Treatment effects by country
quartz(file = 'figs/sim3-te.pdf', type = 'pdf',
       height = 3, width = 4.5)
par(mgp = c(2.2, 1, 0), mar = c(4, 4, 2, 2), cex = .8)
plot(true_te ~ country, data = df_country,
     xlab = 'country', ylab = 'estimated treatment effect', pch = 16, col = "black",
     ylim = c(0, 5),
     xaxt = 'n')
axis(1, at = c(1, 20, 40, 60))
abline(h = coef(mfit_pool)[2], lty = 2)
abline(h = coef(mfit_fe)[2], col = 'tomato')
points(te_msep ~ country, data = df_country, col = 'royalblue', pch = 16)
points(te_mpart ~ country, data = df_country)
text(10, 4.8, 'small n (10)')
text(30, 4.8, 'medium n (50)')
text(50, 4.8, 'large n (100)')
abline(v = c(20.5, 40.5), lwd = .5)
dev.off()


## Figure
# Before matching
intercept_fe <- coef(fit_fe)[1]
fixef <- c(0, coef(fit_fe)[7:65])
intercept_fe <- ifelse(is.na(fixef), intercept_fe, intercept_fe + fixef)
intercept_sep <- coef(fit_sep)[1]
fixef <- c(0, coef(fit_sep)[3:61])
intercept_sep <- ifelse(is.na(fixef), intercept_sep, intercept_sep + fixef)
intercept_part <- coef(fit_part)[1:60]

# one figure for each country
for (i in 1:n_countries) {
    df_one <- filter(df, country == i)
    filename <- paste0('figs/sim2-country-', i, '.pdf')
    quartz(file = filename, type = 'pdf',
           height = 3, width = 4.5)
    par(mgp = c(2.2, 1, 0), mar = c(4, 4, 2, 2), cex = .8)
    plot(outcome ~ treatment, data = df_one, pch = 16, xaxt = 'n',
         main = paste('country', i, sep = ' '))
    axis(1, at = c(0, 1), labels = c('control', 'treatment'))
    abline(a = coef(fit_pool)[1] + 
               coef(fit_pool)[3]*mean(df_one$x1) +
               coef(fit_pool)[4]*mean(df_one$x2) +
               coef(fit_pool)[5]*mean(df_one$z1) +
               coef(fit_pool)[6]*mean(df_one$z2),
           b = coef(fit_pool)[2])
    abline(a = intercept_fe[i] + 
               coef(fit_fe)[3]*mean(df_one$x1) +
               coef(fit_fe)[4]*mean(df_one$x2) +
               coef(fit_fe)[5]*mean(df_one$z1) +
               coef(fit_fe)[6]*mean(df_one$z2),
           b = df_country$te_fe[i],
           col = 'tomato')
    if (!is.na(df_country$te_sep[i])) {
        abline(a = intercept_sep[i] + 
                   coef(fit_sep)[62]*mean(df_one$x1) +
                   coef(fit_sep)[63]*mean(df_one$x2),
               b = df_country$te_sep[i],
               col = 'royalblue')
    }
    abline(a = intercept_part[i] +
               coef(fit_part)[125]*mean(df_one$x1) +
               coef(fit_part)[127]*mean(df_one$x2) +
               coef(fit_part)[126]*mean(df_one$z1) +
               coef(fit_part)[128]*mean(df_one$z2),
            b = df_country$te_part[i],
            lty = 2)
    dev.off()
}          


## After matching
intercept_mfe <- coef(mfit_fe)[1]
fixef <- c(0, coef(mfit_fe)[5:63])
intercept_mfe <- ifelse(is.na(fixef), intercept_mfe, intercept_mfe + fixef)
intercept_msep <- coef(mfit_sep)[1]
fixef <- c(0, coef(mfit_sep)[3:61])
intercept_msep <- ifelse(is.na(fixef), intercept_msep, intercept_msep + fixef)
intercept_mpart <- coef(mfit_part)[1:60]

## One figure for each country
for (i in 1:n_countries) {
    df_one <- filter(mdf, country == i)
    filename <- paste0('figs/sim2-country-m', i, '.pdf')
    quartz(file = filename, type = 'pdf',
           height = 3, width = 4.5)
    par(mgp = c(2.2, 1, 0), mar = c(4, 4, 2, 2), cex = .8)
    plot(outcome ~ treatment, data = df_one, pch = 16, xaxt = 'n',
         main = paste('country', i, 'after matching', sep = ' '))
    axis(1, at = c(0, 1), labels = c('control', 'treatment'))
    abline(a = coef(mfit_pool)[1] + 
               coef(mfit_pool)[3]*mean(df_one$x1) +
               coef(mfit_pool)[4]*mean(df_one$z1),
           b = coef(mfit_pool)[2])
    abline(a = intercept_mfe[i] + 
               coef(mfit_fe)[3]*mean(df_one$x1) +
               coef(mfit_fe)[4]*mean(df_one$z1),
           b = df_country$te_mfe[i],
           col = 'tomato')
    if (!is.na(df_country$te_msep[i])) {
        abline(a = intercept_msep[i] + 
                   coef(mfit_sep)[62]*mean(df_one$x1),
               b = df_country$te_msep[i],
               col = 'royalblue')
    }
    abline(a = intercept_mpart[i] +
               coef(mfit_part)[125]*mean(df_one$x1) +
               coef(mfit_part)[126]*mean(df_one$z1),
           b = df_country$te_mpart[i],
           lty = 2)
    dev.off()
}          
