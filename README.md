# Instructions

This is an implementation of the parameter estimation methods for linear preferential attachment model described in _Fitting the linear preferential attachment model_ (Wan et al.) [1]. The R-scripts `netfnsMLE.R` and `netfnsSnap.R` correspond to estimation algorithms for the full network (MLE) and a snapshot of the network, where the details can be found in Section 3 and 4 of Wan et al. repectively.

### Parameter estimation from the full network (MLE)

This algorithm provide the MLE for the parameters of a 5-scenario linear preferential attachment model given an observed network with full history and an designated time interval for the analysis.

```r
source(‘netfnsMLE.R’)
```

##### Usage

```r
netEstMLE(network,start_time,end_time)
```

##### Input

**network**

The representation of a network with _n_ edges. An _n_ by 3 matrix with each row corresponding to the start node, the end node, and the timestamp of each edge of the network.

**start_time**

The starting point of the time interval to be analyzed. Default: the beginning time of the observations.

**end_time**

The ending point of the time interval to be analyzed. Default: the ending time of the observations.

##### Output

The parameter estimates for \alpha,\beta,\gamma,\xi,\rho,\delta_in,\delta_out.


### Parameter estimation from a network snapshot

This algorithm provide the snapshot estimate for the parameters of a 3-scenario linear preferential attachment model given an observed snapshot of a network.

```r
source(‘netfnsSnap.R’)
```

##### Usage

```r
netEstSnap(network)
```

##### Input

**network**

The representation of a network snapshot with _n_ edges. An _n_ by 2 matrix with each row corresponding to the start node and the end node.


##### Output

The parameter estimates for \alpha,\beta,\gamma,\delta_in,\delta_out.


### Demonstrations

A demonstrative R script and R dataset can be found in `netExample.R` and `netExample.RData`.

### References
[1] https://arxiv.org/abs/1703.03095

