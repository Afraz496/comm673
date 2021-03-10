clc
clear
close all

%The following script implements the Shooting Method in Brunnermeir and
%Sannikov (2016)

%VARIABLE DEFINITIONS
%a: productivity parameter
%iota_t: reinvestment rate per unit of capital
%k_t: physical capital
%fi(iota_t): investment function with adjustment costs
%fi(0)= 0
%K_t: Aggregate amount of capital
%q_t: Price of Capital
%q_t*K_t: Aggregate Net Worth in the economy
%N_t: Net worth of experts
%q_t*K_t - N_t: Net worth of household
%n_t = N_t/(q_t*K_t) in [0,1]: Expert's Wealth share

%Parameter Values
rho = 0.06;
r = 0.05;
a = 0.11;
asmall = 0.05;
delta = 0.03;
chismall = 1;



%Set Terminal Conditions
f = @(iota)((asmall - iota)/(r - fi + deltasmall));
g = @(x)-f(x);

%this line searches for the maximum value
q0 = fminsearch(g,[0 0]);
theta0 = 1;
dtheta0 = -10e-10;
qL = 0;
qH = 10e15;

for i = 1:50
    dq0 = (qL+qH)/2;
    
end