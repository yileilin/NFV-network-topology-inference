# NFV-network-topology-inference
Citation:
Yilei Lin, Ting He, Shiqiang Wang, Kevin S. Chan, and Stephen Pasteris, "Looking Glass of NFV: Inferring the Structure and State of NFV Network from External Observations", IEEE/ACM Transactions on Networking, accepted March 2020.

This work is aim to infer NFV network topology with both single source and multiple sources.
The code is done by matlab and we use Cplex to deal with ILP to reduce running time.

We generate our ground truth using generate_inputs(_multi).m from three dataset (AS1755.mat, AS3967.mat and AS6461.mat).

We save our inferred graph using graph_saver(_multi).m. The rusults include ground truth and graphs inferred by RNJ(REA), CE and SAP.
For simplicity, we prepare part of the matrix used in ILP using SAP_prepare.m.
