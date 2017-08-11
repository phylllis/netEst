# netExample.RData contains a network generated from the linear preferential attachment data with parameter values
# alpha=.2, beta=.5, gamma=.3, xi=rho=0, din=1, dout=2, n=50000
# starting with a self-loop 0->0
# columns:
#   - from: the node the edge pointed from
#   - to: the node the edge pointed to
#   - time: the order of the creation of the edge, served as the timestamp
load("netExample.RData")

# MLE parameter estimates from the full data with the timestamps
source('netfnsMLE.R')
netEstMLE(network)

# MLE parameter estimates from the first half of the data (edge 1 to edge 25000) with the timestamps
source('netfnsMLE.R')
netEstMLE(network,2,25000)

# MLE parameter estimates from the second half of the data (edge 25001 to edge 50000) with the timestamps
source('netfnsMLE.R')
netEstMLE(network,25001,50000)

# Snapshot parameter estimates from the snapshot without the timestamps
source('netfnsSnap.R')
netEstSnap(network[,1:2])

# The search range for d_in and d_out are defaulted as (0.01,100), but can be adjusted through setting din_l, din_u, dout_l, dout_u.
