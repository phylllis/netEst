# The snapshot estimating function
# input:  
#     - network with time stamps
#         - columns: start node, end node
#     - din_l,din_u,dout_l,dout_u: the lower and upper bounds for the numerical search for din and dout
#         - default: lower - 0.1, upper - 100
# output:  parameter estimates: din,dout,alpha,beta,gamma
#         - note that here alpha/beta/gamma are normalized to sum up to 1

netEstSnap <- function(network,din_l=.01,din_u=100,dout_l=.01,dout_u=100){
  colnames(network) <- c('start','end')
  deg_in <- tabulate(network$end+1)
  deg_out <- tabulate(network$start+1)
  
  Ni <- tabulate(deg_in+1) # in-degree distribution N_i^{in}(n_0)
  Ni_over <- c(0,cumsum(rev(Ni))[-length(Ni)]) # N_{>i}^{in}(n_0)
  Ni_ind <- max(deg_in):0
  Nj <- tabulate(deg_out+1) # out-degree distribution N_j^{in}(n_0)
  Nj_over <- c(0,cumsum(rev(Nj))[-length(Nj)]) # N_{>j}^{out}(n_0)
  Nj_ind <- max(deg_out):0
  
  n <- nrow(network)
  N <- length(deg_in)
  p0_in <- Ni[1]/n
  p0_out <- Nj[1]/n  
  
  beta_hat <- 1-N/n

  f1_in <- function(eta){
    1/sapply(eta, function(eta){
      sum(Ni_ind*Ni_over/n/(Ni_ind+eta*(1-Ni_ind*(1-beta_hat))))
    })
  }
  f2_in <- function(eta){
    (1-p0_in*eta)/(beta_hat+p0_in)
  }

  g <- function(eta){
    f1_in(eta)-f2_in(eta)
  }
  eta_in <- uniroot(g, c(0.1, 1/(1-beta_hat)-0.1))$root
  d1_hat <- eta_in/(1-eta_in*(1-beta_hat))
  alpha_hat <- (p0_in+beta_hat)/(1-p0_in*eta_in) - beta_hat

  f1_out <- function(eta){
    1/sapply(eta, function(eta){
      sum(Nj_ind*Nj_over/n/(Nj_ind+eta*(1-Nj_ind*(1-beta_hat))))
    })
  }
  f2_out <- function(eta){
    (1-p0_out*eta)/(beta_hat+p0_out)
  }

  g <- function(eta){
    f1_out(eta)-f2_out(eta)
  }
  eta_out <- uniroot(g, c(0.1, 1/(1-beta_hat)-0.1))$root
  d2_hat <- eta_out/(1-eta_out*(1-beta_hat))
  gamma_hat <- (p0_out+beta_hat)/(1-p0_out*eta_out) - beta_hat

  ag_hat <- alpha_hat+gamma_hat

  # use alpha
  gamma_hat <- gamma_hat/ag_hat*(1-beta_hat)
  alpha_hat <- alpha_hat/ag_hat*(1-beta_hat)

  d2_reest <- function(d2){
    sapply(d2,function(d2){
      sum(Nj_over/(Nj_ind+d2))/n - alpha_hat/d2 -
        (1-alpha_hat)*(1-beta_hat)/(1+(1-beta_hat)*d2)
    })
  }
  d2_hat <- uniroot(d2_reest,c(dout_l,dout_u))$root

  d1_reest <- function(d1){
    sapply(d1,function(d1){
      sum(Ni_over/(Ni_ind+d1))/n - gamma_hat/d1 -
        (1-gamma_hat)*(1-beta_hat)/(1+(1-beta_hat)*d1)
    })
  }
  d1_hat <- uniroot(d1_reest,c(din_l,din_u))$root
  return(list(d_in=d1_hat,d_out=d2_hat,alpha=alpha_hat,beta=beta_hat,gamma=gamma_hat))
}

