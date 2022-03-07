# mixed-traffic-amod-route-rebalance

Solve different problems regarding social and user centring routing in transportation networks while considering:

1. *mixed traffic*: interaction of collaborative (via an Autonomous Mobility-on-Demand system) and user-centric vehicles;
2. *intermodality*: access to different modes of transporation;
3. *rebalancing*: impose rebalaning/empty driving for the AMoD system. 

The specific approach to solve the joint routing and rebalancing problem is detailed in our paper [Wollenstein-Betech et al., 2021](https://ieeexplore.ieee.org/abstract/document/9541261). 
The network and demand data is provided as in [TNTP](https://github.com/bstabler/TransportationNetworks) data format.

* To run use: `python3 -m experiments.{name of the file without extension}`

* Requirements: gurobipy, networkx, scipy, numpy, pwlf
