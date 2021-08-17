# Data-Driven-ILC
Data-Driven ILC with Robustness to Stochastic Parametric Errors

The files implement the methods developed in the paper “Robust Data Driven Iterative Learning Control for Linear-Time-Invariant and Hammerstein-Wiener Systems”, authored by Jianfei Dong, which has been recently accepted for publication in IEEE Transactions on Cybernetics.

The main file of the demo program is “demo_rddilc_hamwien_MC.m”, which calls the other (pre-complied for fast implementation) functions in the same folder. Please note that these are Monte Carlo simulations for a few hundred times. Running it will take some time. The simulations are based on randomly generated stochastic disturbances. Therefore, the results may not be exactly the same with those written in the paper, but shall be close enough.

All the codes are developed in MATLAB R2014b, and shall be put into the same folder.

For suggestions or requests, please contact with the author.

Copyright: Jianfei Dong (C), on Apr. 26, 2021.
