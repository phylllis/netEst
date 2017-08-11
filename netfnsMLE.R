# Step 1: label the nodes according to the order they are created.
# input: network with edges in time order
#   columns: start node, end node, (timestamp)
# output: network with re-labeled nodes according to the order of their creation, starting from 0
nodeLabel <- function(network){
  in_nodes <- network[,2]
  out_nodes <- network[,1]
  tmp <- cbind(out_nodes,in_nodes)
  tmp1 <- matrix(as.numeric(factor(t(tmp),levels=unique(as.vector(t(tmp))))),byrow=T,ncol=2)
  in_nodes <- tmp1[,2]
  out_nodes <- tmp1[,1]
  return(cbind(out_nodes,in_nodes))
}

# Step 2: infer the scenario under which each edge is created
# input: network with ordered nodes
#   columns: start node, end node, (timestamp)
# output: network with two new columns indicating the scenario of the creation of each edge and the number of nodes at the time
#   columns: start node, end node, scenario, # of nodes at the time
nodeEvo <- function(network){
  n_edge <- nrow(network)
  in_nodes <- network[,2]
  out_nodes <- network[,1]
  maxnode <- apply(cbind(cummax(in_nodes),cummax(out_nodes)),1,max)
  t <- 1:n_edge
  evolution <- c(1,diff(maxnode))
  evolution[evolution!=0] <- 1
  evolution[cummax(in_nodes)>cummax(out_nodes)&evolution!=0] <- 2
  evolution[cummax(in_nodes)>c(0,maxnode[-n_edge])&
              cummax(out_nodes)>c(0,maxnode[-n_edge])] <- 3
  evolution[cummax(in_nodes)>c(0,maxnode[-n_edge])&
              cummax(out_nodes)==cummax(in_nodes)] <- 4
  # evolution: 0-beta; 1-alpha; 2-gamma; 3-new-new
  evo <- rep(NA,n_edge)
  evo[evolution==1] <- 1
  evo[evolution==0] <- 2
  evo[evolution==2] <- 3
  evo[evolution==3] <- 4
  evo[evolution==4] <- 5
  Nt_ct <- rep(0,n_edge)
  Nt_ct[evolution==1] <- 1
  Nt_ct[evolution==0] <- 0
  Nt_ct[evolution==2] <- 1
  Nt_ct[evolution==3] <- 2
  Nt_ct[evolution==4] <- 1
  Nt <- cumsum(Nt_ct)
  return(cbind(out_nodes,in_nodes,evo,Nt))
}

# Step 3: estimating the parameters during the designated time interval
# input: 
#     - network with evolution and number of nodes from Step 2
#         - columns: start node, end node, scenario, # of nodes at the time
#     - start_edge, end_edge: the start and end edge index between which the model is to be estimated (start_edge>=2)
#     - din_l,din_u,dout_l,dout_u: the lower and upper bounds for the numerical search for din and dout
# output: parameter estimates: din,dout,alpha,beta,gamma,xi,rho
parEst <- function(network,start_edge,end_edge,din_l,din_u,dout_l,dout_u){
  in_nodes_start <- network[1:(start_edge-1),2]
  out_nodes_start <- network[1:(start_edge-1),1]
  in_nodes_end <- network[1:end_edge,2]
  out_nodes_end <- network[1:end_edge,1]
  deg_in_start <- tabulate(in_nodes_start)
  deg_out_start <- tabulate(out_nodes_start)
  deg_in_end <- tabulate(in_nodes_end)
  deg_out_end <- tabulate(out_nodes_end)
  ######
  Ni_start <- tabulate(deg_in_start+1) # in-degree distribution N_i^{in}(n_0)
  Ni_over_start <- c(0,cumsum(rev(Ni_start))[-length(Ni_start)]) # N_{>i}^{in}(n_0)
  Ni_end <- tabulate(deg_in_end+1) # in-degree distribution N_i^{in}(n_0)
  Ni_over_end <- c(0,cumsum(rev(Ni_end))[-length(Ni_end)]) # N_{>i}^{in}(n_0)
  Ni_over_start <- c(rep(0,length(Ni_over_end)-length(Ni_over_start)),Ni_over_start)
  Mi_over <- Ni_over_end - Ni_over_start
  i_ind <- max(deg_in_end):0
  ######
  Nj_start <- tabulate(deg_out_start+1) # out-degree distribution N_j^{out}(n_0)
  Nj_over_start <- c(0,cumsum(rev(Nj_start))[-length(Nj_start)]) # N_{>j}^{out}(n_0)
  Nj_end <- tabulate(deg_out_end+1) # out-degree distribution N_j^{out}(n_0)
  Nj_over_end <- c(0,cumsum(rev(Nj_end))[-length(Nj_end)]) # N_{>j}^{out}(n_0)
  Nj_over_start <- c(rep(0,length(Nj_over_end)-length(Nj_over_start)),Nj_over_start)
  Mj_over <- Nj_over_end - Nj_over_start
  j_ind <- max(deg_out_end):0
  ######
  Nt <- network$Nt[(start_edge):end_edge]
  evolution <- network$evo[(start_edge):end_edge]
  Nt_in <- Nt[evolution %in% c(1,2)]
  t_in <- ((start_edge):end_edge)[evolution %in% c(1,2)]
  Nt_out <- Nt[evolution %in% c(2,3)]
  t_out <- ((start_edge):end_edge)[evolution %in% c(2,3)]
  ######
  d2_left <- function(d2){
    sapply(d2,function(d2){sum(Mj_over/(j_ind+d2))}) - sum(evolution %in% c(1,4,5))/d2
  }
  d2_right <- function(d2){
    sapply(d2,function(d2){sum(Nt_out/(t_out + d2*Nt_out))})
  }
  f2 <- function(d2){
    d2_left(d2) - d2_right(d2)
  }
  d2_hat <- uniroot(f2,c(dout_l,dout_u))$root
  ######
  d1_left <- function(d1){
    sapply(d1,function(d1){sum(Mi_over/(i_ind+d1))}) - sum(evolution %in% c(3,4,5))/d1
  }
  d1_right <- function(d1){
    sapply(d1,function(d1){sum(Nt_in/(t_in + d1*Nt_in))})
  }
  f1 <- function(d1){
    d1_left(d1) - d1_right(d1)
  }
  d1_hat <- uniroot(f1,c(din_l,din_u))$root
  alpha_hat <- sum(evolution==1)/(end_edge-start_edge+1)
  beta_hat <- sum(evolution==2)/(end_edge-start_edge+1)
  gamma_hat <- sum(evolution==3)/(end_edge-start_edge+1)
  xi_hat <- sum(evolution==4)/(end_edge-start_edge+1)
  rho_hat <- sum(evolution==5)/(end_edge-start_edge+1)
  return(list(d_in=d1_hat,d_out=d2_hat,alpha=alpha_hat,beta=beta_hat,gamma=gamma_hat,
              xi=xi_hat,rho=rho_hat))
}

# The MLE estimating function
# input: 
#     - network with time stamps
#         - columns: start node, end node, timestamp
#     - start and end time stamps of the time interval to be estimated
#         - default: the whole period observed
#     - din_l,din_u,dout_l,dout_u: the lower and upper bounds for the numerical search for din and dout
#         - default: lower - 0.1, upper - 100
# output:  parameter estimates: din,dout,alpha,beta,gamma,xi,rho
netEstMLE <- function(network,start_time=min(network[,3]),end_time=max(network[,3]),
                      din_l=.01,din_u=100,dout_l=.01,dout_u=100){
  colnames(network) <- c('start','end','time')
  net <- network[order(network$time),]
  net1 <- nodeLabel(net)
  remove(net)
  net2 <- nodeEvo(net1)
  remove(net1)
  net2 <- data.frame(net2)
  start_edge <- sum(network$time<=start_time)+1
  end_edge <- sum(network$time<end_time)+1
  par_est <- parEst(net2,start_edge,end_edge,din_l,din_u,dout_l,dout_u)
  return(par_est)
}

