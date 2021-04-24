function [Phi iota] = investment(q)
% [Phi iota] = investment(q)
% finds the optimal investment rate iota and rate of capital creation Phi
% given the price q

kappa = 10;

theta = 10;  % adjustment cost parameter
iota = (q - 1)/(2*kappa - 1);  Phi = (1/kappa)*(sqrt(1+2*kappa*iota) - 1);
