require(glmmTMB)

simulate_data <- function(model,
                          nsim=1) {
  if(nsim>1) {
    sim_data <- matrix(NA, ncol=nsim, nrow=nrow(model$frame))
  } else {
    sim_data <- rep(NA, nrow(model$frame))
  }

  if(model$modelInfo$family$family == "poisson") {
    mu <- predict(model, type="response", re.form = NULL)

    if(is.matrix(sim_data)) {
      for(i in 1:nsim) {
        sim_data[,i] <- rpois(length(mu), lambda=mu)
      }
    } else {
      sim_data <- rpois(length(mu), lambda=mu)
    }
  } else if(model$modelInfo$family$family == "nbinom2") {
    mu <- predict(model, type="response", re.form = NULL)
    size <- sigma(model)

    if(is.matrix(sim_data)) {
      for(i in 1:nsim) {
        sim_data[,i] <- rnbinom(length(mu), mu=mu, size=size)
      }
    } else {
      sim_data <- rnbinom(length(mu), mu=mu, size=size)
    }
  }

  return(list(mu=mu, sim=sim_data))
}

