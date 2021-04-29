# COMM 673: Advanced Topics in Theoretical Asset Pricing at UBC

This Repository holds all the code written for the paper replication for the 25% grade requirement at UBC.

The paper we chose to replicate is Brunnermeir and Sannikov (2016) which can be found [here](https://scholar.princeton.edu/sites/default/files/markus/files/_macrohandbook_brusan.pdf)

This repository consists of 3 code files:

1. `simplemodel.m`: Replicates the plots on page 16 from Section 2 of the paper
2. `shootingplots.m`: Implements the Shooting Method which is under Section 3.4 Method 1 (page 27). It also plots the change of sigma.
3. `plotiterative.m`: Implements the Iterative Method which is under Section 3.5 Method 2 (page 33) on some key parameters of interest.

## Simple Model

**Note**: We are not including this as our submission, we believe this wasn't really a replication but we have left the code for it.

Quite an easy replication. Parameters of interest are:

- $a = 0.1$
- $\rho = 5%$
- $\sigma = 0.1$
- $\fi = \frac{log(\kappa\iota + 1)}{\kappa}$
- $\kappa = 10$

The following Values were extrapolated from different parts of the paper:

- $\iota = \frac{(q-1)}{\kappa}$ -> Implies $\iota = 0.04$
- $\delta = 0.0336$: Note that this is exactly equal to $\fi(0.04)$

The full MATLAB script prints the plots in their entirety on page 16.

## Shooting Method

Our Shooting Method borrows the sample code from Brunnermeir and Sannikov (2014) and modifies their investment function and some other variables to generate the output in Brunnermeir and Sannikov (2016). To run our file navigate to `112732 > data`. The script of interest is `shootingplots.m` and if the MATLAB path is correct, it should run fine.

## Iterative Method

We would like to thank Professor Lorenzo for sharing the Iterative Method starter code. We modified the files accordingly to generate output for parameter changes and also fixed a bug in the original code which used the productivitiy parameter of experts in both value functions. To run this file in the `main` folder, you must run `plotiterative.m` which will generate a MATLAB figures output of the parameters changed under the Value function iteration section of the paper. It works well for the first 2 parameters which require change (`sigma` and `chibar`) but it runs into convergence issues with the productivity parameter and the risk aversion.
