testDistribution <- function(p, x, distribution = "normal", giveDistributions = F){
  if (giveDistributions == T){
    NLL <- c("normal", "gamma", "beta", "negative binomial")
    distribution = "none"
  }
  if (str_to_lower(distribution) == "normal"){
    mu <- p[1]
    sigma <- p[2]
    NLL <- -sum(dnorm(x, mean = mu, sd = sigma, log = T))
  }
  if (str_to_lower(distribution) == "gamma"){
    shape <- p[1]
    rate <- p[2]
    NLL <- -sum(dgamma(x, shape = shape, rate = rate, log = T))
  }
  if (str_to_lower(distribution) == "beta"){
    #The beta distribution only ranges from [0;1] and thus it is
    #exclusively relevant to the normalized wind power statistic
    #and not the fitting of the other two parameters.
    shape1 <- p[1]
    shape2 <- p[2]
    NLL <- -sum(dbeta(x, shape1 = shape1, shape2 = shape2, log = T))
  }
  if (str_to_lower(distribution) == "exponential"){
    lambda = p
    NLL <- -sum(dexp(x, rate = lambda, log = T))
  }
  
  if (str_to_lower(distribution) == "weibull"){
    shape = p[1]
    scale = p[2]
    NLL <- -sum(dweibull(x, shape = shape, scale = scale, log = T))
  }
  if (str_to_lower(distribution) == "negative binomial"){
    #Virker ikke endnu
    alpha <- p[1] #target number of succesfull trials
    probs <- p[2] #probability of succes in each trial
    NLL <- -sum(dnbinom(x = alpha, size = x, prob = probs, log = T))
  }
  if (str_to_lower(distribution) == "binomial"){
    p <- p
    n <- x[2]
    k <- x[1]
    NLL <- -sum(dbinom(x = k, size = n, prob = p, log = T))
  }
  if (str_to_lower(distribution) == "poisson"){
    lambda = p
    NLL <- -sum(dpois(x = x, lambda = lambda, log = F))
  }
  return(NLL)
}
