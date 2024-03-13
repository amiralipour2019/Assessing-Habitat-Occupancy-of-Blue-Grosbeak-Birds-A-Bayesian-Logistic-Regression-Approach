# Assuming we have a dataset `data` with fields `Y` (binary outcomes), `X` (covariates), and `ns` (number of fields)
bg_data=read.csv("Bayesian_logistic_regression.csv")[,1:5]
head(bg_data)

# Set up the covariate (field size) and outcome matrix
X <- bg_data$field.size
Y <- bg_data[, 2:4]  # Assuming the binary outcomes are in columns 2 to 4

set.seed(56) # For reproducibility
################################################################### Initialization:
#Set up the covariate (field size) and outcome matrix
X=bg_data$field.size
Y=bg_data[,2:4]

set.seed(56)

# Initial parameter values
alpha_current=0
beta_current=0
theta_current=0.5
c2=1

#Proposal standard deviations for the random walk
sd_alpha=0.1
sd_beta=0.1 
sd_theta=0.05

#Number of iterations for the Metropolis-Hastings algorithm
n_iter=10000
burn_in=5000

#Store samples
samples=matrix(NA,nrow=n_iter-burn_in,ncol=3)
colnames(samples)=c("alpha","beta","theta")
#############################################################
#Function to calculate log likelihood
calculate_log_likelihood=function(alpha,beta,theta,X,Y){
  
  log_likelihood=0
  
  ns=length(X) 
  nt=ncol(Y)  
  
  for (s in 1:ns){
    for (t in 1:nt){
      
      y_st=Y[s,t]
      x_s=X[s] 
      z_s_prob=1/(1+exp(-(alpha+x_s*beta)))  
      p_y_given_z_theta=(theta*z_s_prob)^y_st*(1-theta*z_s_prob)^(1-y_st)
      log_likelihood=log_likelihood+log(p_y_given_z_theta)
    }
  }
  return(log_likelihood)
}

# Generate samples
for (i in 1:n_iter){
  
  #Propose new values
  alpha_proposed=rnorm(1,alpha_current,sd_alpha)
  beta_proposed=rnorm(1,beta_current,sd_beta) 
  theta_proposed=runif(1, max(0,theta_current-sd_theta),min(1,theta_current+sd_theta))
  
  # Calculate log likelihood for current and proposed
  log_likelihood_current=calculate_log_likelihood(alpha_current,beta_current,theta_current,X,Y)
  log_likelihood_proposed=calculate_log_likelihood(alpha_proposed,beta_proposed,theta_proposed,X,Y)
  
  # Calculate log prior for current and proposed
  log_prior_current= -0.5/c2*(alpha_current^2+beta_current^2)+log(dunif(theta_current,0,1))
  log_prior_proposed= -0.5/c2*(alpha_proposed^2+beta_proposed^2)+log(dunif(theta_proposed,0,1))
  
  # Calculate log posterior for current and proposed
  log_posterior_current=log_likelihood_current+log_prior_current
  log_posterior_proposed=log_likelihood_proposed+log_prior_proposed
  
  # Acceptance ratio
  r=exp(log_posterior_proposed-log_posterior_current)
  
  # Accept or reject the proposed values
  if (runif(1)<r){
    alpha_current=alpha_proposed
    beta_current=beta_proposed
    theta_current=theta_proposed
  }
  
  if (i>burn_in){
    samples[i-burn_in,]=c(alpha_current,beta_current,theta_current)
  }
}


head(samples)
###############################################################################
# Trace Plots
par(mfrow=c(3,1))  
# Plotting alpha

plot(samples[, "alpha"],type= "l",main="Trace Plot of Alpha", xlab="Iteration",ylab="Alpha")

# Plotting beta
plot(samples[, "beta"],type="l",main="Trace Plot of Beta", xlab="Iteration",ylab="Beta")

# Plotting theta
plot(samples[, "theta"],type="l", main="Trace Plot of Theta", xlab="Iteration", ylab ="Theta")
#########################
# Beta distribution
hist(samples[,"beta"], breaks=40, main="Posterior Distribution of Beta",xlab ="Beta")
abline(v =0,col = "red")

# Derive the intervals
beta_samples=samples[,"beta"]
beta_cred_int=quantile(beta_samples,probs = c(0.025, 0.975))
print(paste("95% Credible Interval: [", beta_cred_int[1], ", ",beta_cred_int[2], "]",sep=""))

# Let define the function to compute estimated occupancy probability from model parameters

calculate_occupancy_probability= function(alpha, beta,field_size){
  
  1/(1 + exp(-(alpha+beta*field_size)))
}

# Set threshold 
threshold=0.2


# Store the values
occupancy_probabilities=matrix(NA,nrow=length(X),ncol=n_iter-burn_in)

for (iter in 1:n_iter-burn_in) {
  for (s in 1:length(X)) {
    occupancy_probabilities[s, iter]=calculate_occupancy_probability(
      alpha=samples[iter, "alpha"],
      beta=samples[iter, "beta"],
      field_size=X[s]
    )
  }
}


estimated_occupancies=occupancy_probabilities>threshold


posterior_probability_all_occupied=mean(apply(estimated_occupancies, 2,all))

#output

print(posterior_probability_all_occupied)

# Assess the sensitivity of c
library(rstanarm)
c2_values=c(0.1, 1, 10, 100) 
results=list()

for (c2 in c2_values) {
  
  model_fit=stan_glm(
    
    bg_data[, 2] ~ bg_data$field.size,
    
    family=binomial(link = "logit"),
    
    prior=normal(0, sqrt(c2)), 
    
    prior_intercept=normal(0, sqrt(c2)),
    
    iter=2000,
    
    seed=56
  )
  
 
  
  results[[as.character(c2)]]=summary(model_fit)
 
}
