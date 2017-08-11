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
  deg_in<- tabulate(network$end)
  deg_out <- tabulate(network$start)
  
  Ni <- tabulate(deg_in+1) # in-degree distribution N_i^{in}(n_0)
  Ni_over <- c(0,cumsum(rev(Ni))[-length(Ni)]) # N_{>i}^{in}(n_0)
  i_ind <- max(deg_in):0
  Nj <- tabulate(deg_out+1) # out-degree distribution N_j^{in}(n_0)
  Nj_over <- c(0,cumsum(rev(Nj))[-length(Nj)]) # N_{>j}^{out}(n_0)
  j_ind <- max(deg_out):0
  
  n <- nrow(network)
  beta_hat <- 1-length(unique(c(network$start,network$end)))/n
  d1_left <- function(d1){
    sapply(d1,function(d1){
      sum(Ni_over*i_ind/(i_ind+d1))*(1+d1*(1-beta_hat))/n
    })
  }
  d1_right <- function(d1){
    (Ni[1]/n+beta_hat)/(1-Ni[1]/n*d1/(1+d1*(1-beta_hat)))
  }
  f1 <- function(d1){
    d1_left(d1) - d1_right(d1)
  }
  d1_hat <- uniroot(f1,c(din_l,din_u))$root
  alpha_hat <- d1_right(d1_hat)-beta_hat
  
  d2_left <- function(d2){
    sapply(d2,function(d2){
      sum(Nj_over*j_ind/(j_ind+d2))*(1+d2*(1-beta_hat))/n
    })
  }
  d2_right <- function(d2){
    (Nj[1]/n+beta_hat)/(1-Nj[1]/n*d2/(1+d2*(1-beta_hat)))
  }
  f2 <- function(d2){
    d2_left(d2) - d2_right(d2)
  }
  d2_hat <- uniroot(f2,c(dout_l,dout_u))$root
  gamma_hat <- d2_right(d2_hat)-beta_hat
  
  
  ag_hat <- alpha_hat+gamma_hat
  # use alpha
  gamma_hat <- gamma_hat/ag_hat*(1-beta_hat)
  alpha_hat <- alpha_hat/ag_hat*(1-beta_hat)
  
  f2_reest <- function(d2){
    sapply(d2,function(d2){
      sum(Nj_over/(j_ind+d2))/n - alpha_hat/d2 -
        (1-alpha_hat)*(1-beta_hat)/(1+(1-beta_hat)*d2)
    })
  }
  d2_hat <- uniroot(f2_reest,c(dout_l,dout_u))$root
  
  f1_reest <- function(d1){
    sapply(d1,function(d1){
      sum(Ni_over/(i_ind+d1))/n - gamma_hat/d1 -
        (1-gamma_hat)*(1-beta_hat)/(1+(1-beta_hat)*d1)
    })
  }
  d1_hat <- uniroot(f1_reest,c(din_l,din_u))$root
  return(list(d_in=d1_hat,d_out=d2_hat,alpha=alpha_hat,beta=beta_hat,gamma=gamma_hat))
}

