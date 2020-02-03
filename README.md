# mixed-traffic-amod-route-rebalance

Solve different problems regarding social and selfish routing in transportation networks. 
Also, uses the CARS model (see [Salazar et al., 2019](http://asl.stanford.edu/wp-content/papercite-data/pdf/Salazar.Tsao.ea.ECC19.pdf)) to solve the problem of jointly routing and rebalancing AMoD services. 
The network and demand data is provided as in [TNTP](https://github.com/bstabler/TransportationNetworks) data format.

* To run use: `python3 -m experiments.{name of the file without extension}

* Requirements: gurobipy, networkx, scipy, numpy, pwlf
