# COMM 673: Advanced Topics in Theoretical Asset Pricing at UBC

This Repository holds all the code written for the paper replication for the 25% grade requirement at UBC.

The paper we chose to replicate is Brunnermeir and Sannikov (2016) which can be found [here](https://scholar.princeton.edu/sites/default/files/markus/files/_macrohandbook_brusan.pdf)

This repository consists of 3 code files:

1. `simplemodel.m`: Replicates the plots on page 16 from Section 2 of the paper
2. `shooting.m`: Implements the Shooting Method which is under Section 3.4 Method 1 (page 27). The full algorithm description can be found on the paper.
3. `iterative.m`: Implements the Iterative Method which is under Section 3.5 Method 2 (page 33).

## Simple Model

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

## Iterative Method
