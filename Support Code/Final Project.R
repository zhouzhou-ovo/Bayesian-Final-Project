library(ggplot2)
library(GGally)
library(rms)
library(car)
library(skimr)
library(rstanarm)
library(ResourceSelection)
library(MASS)
library(pROC)
library(mvtnorm)
library(coda)

######## load dataset
White <- read.csv("winequality-white.csv",sep=";")
##### statistical description
### scatterplot matrix
ggpairs(White)

##### statistical description:mean, sd, quantile, quality
skim(White)
table(White$quality)

##### Classification for high quality and low quality
32
subset_white <- White[White$quality >= 6, ]
set.seed(123)
resample_high <- subset_white[sample(1:nrow(subset_white), size = length(which(White$quality < 6)), replace = FALSE), ]
white <- rbind(resample_high,White[White$quality < 6, ])
##### if quality >=5, high quality with Quality=1; otherwise Quality=0
white$Quality <- ifelse(white$quality>=6, 1,0)
white[,1:ncol(white)] <- scale(white[,1:ncol(white)]) # standardize the sample data
white$Quality <- ifelse(white$Quality>0,1,0) # sort high quality and low quality
table(white$Quality)

###### initialize training set and testing set
set.seed(1)
train_index <- sample(1:nrow(white),size = 0.7*nrow(white))
train <- white[train_index,]
test <- white[-train_index,]


######## model fitting
##### WWaFD
model <- glm(Quality~fixed.acidity+ volatile.acidity+ citric.acid+
                      residual.sugar+ chlorides+ free.sulfur.dioxide+
                      total.sulfur.dioxide+ density+ pH+ sulphates+
                      alcohol, data = train, family = binomial)
summary(model)
vif(model) # exist multicollinearity

##### WWaBD
#### our own codes
y <- train$Quality
x <- cbind(rep(1,length(y)),as.matrix(train[1:11])) # Design matrix
x_reduced <- x[,2:12] # Design matrix without first column
p <- 12 # number of regression coefficients
n <- length(y) # number of training samples

pmn.beta <- rep(0,p) # prior mean
psd.beta <- diag(10,p) # prior SD

ACR<-NULL # acceptance ratio
seff <- NULL # effective sample size

delta2 <- 0.025^2 # proposal variance
var.prop <- diag(delta2,nrow = p,ncol = p) # proposal covariance matrix
beta <- rep(0,p) # initialize beta
S <- 1.5*10^4 # number of iteration
BETA <- matrix(0,nrow=S/3,ncol=p) # Thinning
ac <- 0 # initialize acceptance
set.seed(123)

### Metropolis
for(s in 1:S) {
  # propose a new alpha.beta
  beta.star <- t(rmvnorm(1, beta, var.prop))
  
  # log r
  logr <- sum(y*(x%*%beta.star)) - 
    sum(log(1+exp(x%*%beta.star))) - 
    sum(y*(x%*%beta)) + 
    sum(log(1+exp(x%*%beta))) + 
    sum(dmvnorm(as.vector(beta.star),pmn.beta,psd.beta,log=T)) -
    sum(dmvnorm(as.vector(beta),pmn.beta,psd.beta,log=T))
  
  if( log(runif(1)) < logr ) { beta <- beta.star; ac <- ac+1 }
  if(s %% 3 ==0) {BETA[s/3,] <- beta}}

ACR <- c(ACR,ac/S) # compute acceptance ratio

### Traceplot of coefficients
par(mfrow=c(3,2))
for (i in 1:p) {
  traceplot(as.mcmc(BETA[,i]),main = bquote(Traceplot: beta[.(i-1)]))
}

### ACF plot of coefficients
for (i in 1:p) {
  acf(as.mcmc(BETA[,i]),main = bquote(ACF: beta[.(i-1)]))
}

### effective sample size
seff <- effectiveSize(BETA)

######## model selection
#### WWaFD
## stepwise selection
step_model <- step(model)
summary(step_model)

#### WWaBD
## our own codes: Metropolis-Hasting
gamma_prior <- rep(0.5,p-1) # Prior probability for gamma

## Function to compute the log-likelihood
log_likelihood <- function(beta0,beta, gamma, y, X) { 
  theta <- beta0 + X %*% (gamma * beta)
  p <- exp(theta) / (1 + exp(theta))
  sum(y * log(p) + (1 - y) * log(1 - p))
}

## Function to compute the log-prior
log_prior <- function(beta0,beta, gamma) { 
  log_prior_beta0 <- dnorm(beta0, mean = 0,sd = 10, log = TRUE) 
  
  log_prior_beta <- sum(dnorm(beta, mean = 0,sd = 10, log = TRUE) * gamma) 
  
  log_prior_gamma <- sum(dbinom(gamma, size = 1, prob = gamma_prior, log = TRUE))
                                                                                                                                         
  log_prior_beta0 + log_prior_beta + log_prior_gamma 
}

S <- 7000 # number of iteration

## initialize storage for samples
beta0_samples <- numeric(S)
beta_samples <- matrix(0, nrow = S, ncol = p-1) 
gamma_samples <- matrix(0, nrow = S, ncol = p-1) 
accept_beta <- accept_gamma <- 0

## Initial values
beta0 <- 0
beta <- rep(0, p-1) 
gamma <- rep(1, p-1)

## Proposal standard deviations
beta0_sd <- 0.035
beta_sd <- rep(0.035, p-1) 
gamma_prob <- 1 / (p-1)

set.seed(123)

## Run Metropolis-Hastings algorithm
for (s in 1:S) {
  
  # Update gamma 
  gamma_proposed <- gamma
  
  for (j in 1: (p-1)) {
    if (runif(1) < gamma_prob) {
      gamma_proposed[j] <- 1 - gamma[j] 
      }
  }
  
  # Compute log-posterior for current and proposed gamma
  log_post_current <- log_likelihood(beta0, beta, gamma, y, x_reduced) + 
    log_prior(beta0, beta, gamma)
  
  log_post_proposed <- log_likelihood(beta0, beta, gamma_proposed, y, x_reduced) + 
    log_prior(beta0, beta, gamma_proposed)
  
  # Accept/reject gamma
  if (log(runif(1)) < (log_post_proposed - log_post_current)) { 
    gamma <- gamma_proposed
    accept_gamma <- accept_gamma + 1 
  }
  
  # Update beta
  beta0_proposed <- rnorm(1, beta0, beta0_sd) 
  beta_proposed <- beta
  
  for (j in 1:(p-1)) {
    if (gamma[j] == 1) {
        beta_proposed[j] <- rnorm(1, mean = beta[j], sd = beta_sd[j]) 
      }
  }
  
  # Compute log-posterior for current and proposed beta
  log_post_current <- log_likelihood(beta0, beta, gamma, y, x_reduced) + 
    log_prior(beta0, beta, gamma)

  log_post_proposed <- log_likelihood(beta0_proposed, beta_proposed, gamma, y, x_reduced) + 
    log_prior(beta0_proposed, beta_proposed, gamma)
  
  # Accept/reject beta
  if (log(runif(1)) < (log_post_proposed - log_post_current)) { 
    beta0 <- beta0_proposed
    beta <- beta_proposed
    accept_beta <- accept_beta + 1
  }
  
  # Store samples
  beta0_samples[s] <- beta0
  beta_samples[s, ] <- beta
  gamma_samples[s, ] <- gamma
}

## Compute acceptance rates
accept_rate_gamma <- accept_gamma / S 
accept_rate_beta <- accept_beta / S

## sample efficient size
seff_gamma <- effectiveSize(gamma_samples)
seff_beta <- effectiveSize(beta_samples)
seff_beta0 <- effectiveSize(beta0_samples)

## beta*gamma posterior mean
bg_mean <- colMeans(beta_samples*gamma_samples)

## traceplot of beta
par(mfrow = c(3, 2))
for (j in 1:(p-1)) {
  traceplot(as.mcmc(beta_samples[, j]),
            main = bquote(Traceplot: beta[.(j)]))
}

## traceplot of beta*gamma
par(mfrow=c(3,2))
for (j in 1:(p-1)) {
  traceplot(as.mcmc(beta_samples[,j]*gamma_samples[,j]),
            main = bquote(Traceplot: beta[.(j)]*gamma[.(j)]))
}

## acf plot of beta*gamma
par(mfrow = c(3, 2)) 
for (j in 1:(p-1)) {
  acf(as.mcmc(beta_samples[, j] * gamma_samples[, j]), 
      main = bquote(ACF: beta[.(j)] * gamma[.(j)]))
}


## get top 5 most frequently occurring values of gamma
gamma_configurations <- apply(gamma_samples, 1, paste, collapse = "") 
gamma_freq <- sort(table(gamma_configurations), decreasing = TRUE) 
top_5_gamma <- names(gamma_freq)[1:5]
top_5_probs <- gamma_freq[1:5] / S
cat("Top 5 gamma configurations:", top_5_gamma, "\n")
cat("Posterior probabilities of top 5 gamma:",top_5_probs,"\n")

## get the best model by Bayesian Model Selection
top1 <- as.numeric(strsplit(top_5_gamma[1],"")[[1]])
beta_mean <- colMeans(beta_samples)
beta0_mean <- mean(beta0_samples)
beta_selected <- top1 * beta_mean # estimate of coefficients after model selection

## compute p(y=1|gamma,beta,x)
p_y <- exp(beta0_mean + x_reduced %*% beta_selected)
y_prob_fitted <- p_y/(1+p_y)
quality_classes_bayes <- ifelse(y_prob_fitted>0.5,1,0) # classification based on fitted probabilities

## compute confusion matrix based on bayesian analysis
confusion_matrix_bayes <- table(Predicted=quality_classes_bayes,Actual=train$Quality)
print(confusion_matrix_bayes)

## compute accuracy of training set based on bayesian analysis
accuracy_bayes <- mean(quality_classes_bayes==train$Quality)
print(accuracy_bayes)

######## model diagnostics
### Measure of Influence
vif(step_model)
model <- glm(Quality~ volatile.acidity+residual.sugar+free.sulfur.dioxide+pH+sulphates+alcohol, data = train, family = binomial)
summary(model)
vif(model)

### hosmer-lemeshow goodness of fit test
hoslem.test(train$Quality,model$fitted.values,g=10)

### roc curve and auc
pred_probs <- predict(model, type = "response")
roc_curve <- roc(train$Quality, pred_probs)
auc_value <- auc(roc_curve)
par(mfrow=c(1,1))
plot(roc_curve,main = paste("ROC Curve (AUC =", round(auc_value, 2), ")"),
     col = "blue", lwd = 2)

### compute confusion matrix based on frequency
quality_classes_freq <- ifelse(pred_probs>0.5,1,0)
confusion_matrix_freq <- table(Predicted=quality_classes_freq,Actual=train$Quality)
print(confusion_matrix_freq)

### compute accuracy of training set based on frequency
accuracy_freq <- mean(quality_classes_freq==train$Quality)
print(accuracy_freq)

##### Compare prediction of two approaches on testing set
### Frequency
pred_probs_test <- predict(model,newdata = test,type = "response")

# roc curve and auc for testing set
roc_curve <- roc(test$Quality, pred_probs_test)
auc_value <- auc(roc_curve)
par(mfrow=c(1,1))
plot(roc_curve,main = paste("ROC Curve (AUC =", round(auc_value, 2), ")"),
     col = "blue", lwd = 2)

# confusion matrix for testing set based on Frequency
quality_classes_freq_test <- ifelse(pred_probs_test>0.5,1,0)
confusion_matrix_freq_test <- table(Predicted=quality_classes_freq_test,Actual=test$Quality)
print(confusion_matrix_freq_test)

# compute accuracy of testing set
accuracy_freq_test <- mean(quality_classes_freq_test==test$Quality)
print(accuracy_freq_test)

### Bayesian
## compute p(y=1|gamma,beta,x)
x_test <- as.matrix(test[1:11])
p_y_test <- exp(beta0_mean + x_test %*% beta_selected)
y_prob_pred_test <- p_y_test/(1+p_y_test)

## classification based on predicted probabilities
quality_classes_bayes_test <- ifelse(y_prob_pred_test>0.5,1,0) 

## confusion matrix for testing set based on Bayesian
confusion_matrix_bayes_test <- table(Predicted=quality_classes_bayes_test,Actual=test$Quality)
print(confusion_matrix_bayes_test)

## compute accuracy of testing set
accuracy_bayes_test <- mean(quality_classes_bayes_test==test$Quality)
print(accuracy_bayes_test)
