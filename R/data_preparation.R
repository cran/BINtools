
compute_sufficient_statistics <- function(event)
{
  outcome <- event$outcome
  x <- event$control_probits
  y <- event$treatment_probits
  N_0 <- length(x)
  N_1 <- length(y)
  M_0 <- sum(x)
  M_1 <- sum(y)
  V_0 <- sum(x^2)
  V_1 <- sum(y^2)
  C_0 <- sum(outer(x, x, '*')) - V_0
  C_1 <- sum(outer(y, y, '*')) - V_1
  C_01 <- sum(outer(x, y, '*'))
  list(outcome = outcome, N_0 = N_0, N_1 = N_1, M_0 = M_0, M_1 = M_1, V_0 = V_0, V_1 = V_1, C_0 = C_0, C_1 = C_1, C_01 = C_01)
}

summariser <- function(dataset)
{
  suffstat_data <- dplyr::bind_rows(lapply(dataset, compute_sufficient_statistics))
  stan_data <- as.list(suffstat_data)
  stan_data$N <- nrow(suffstat_data)
  list(
    data = list(
      raw = dataset,
      stan_data = stan_data,
      detailed = paste("$\\rho_o =", as.character(attr(dataset, "rho_o")), "$"),
      prefix_reg = "",
      prefix_super = ""
    )
  )
}

data_preparation<-function(Outcomes, Control, Treatment=NULL){

  predictions<-c(unlist(Control),unlist(Treatment))
  if( any(predictions<= 0 | predictions>= 1) ){stop('All predictions should be strictly between 0 and 1.')}
  if(!all(Outcomes %in% 0:1)){stop('All outcomes should be binary.')}

  N=length(Outcomes)

  dd<-list()
  for(n in 1:N){
    dd[[n]] <- list(outcome = Outcomes[[n]],
                    control_probits = qnorm(Control[[n]]),
                    treatment_probits = qnorm(as.numeric(Treatment[[n]]))
                    )
  }
  data = lapply(list(dd), summariser)
  return(data)
}

