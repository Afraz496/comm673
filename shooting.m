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
K = 10;

%In the paper on page 29 they state chi = chismall for all time
%They also state that households demand no riskpremium (Sigma: note the
%capital S, it is the hook symbol)
chi = chismall;
Sigmasmall = 0;

% Modify the Volatility accordingly
sigma = 0.25;
 
%The investment function with adjustment costs
fi = @(iota)((1/K)*(sqrt(1+2*K*iota) - 1));

%Set Terminal Conditions
f = @(iota)((asmall - iota)/(r - fi + deltasmall));
g = @(x)-f(x);

%this line searches for the maximum value over the specified interval of
%iota (we need to constrain iota in the range []
q = fminsearch(g,[0 0]);
theta = 1;
dtheta = -10e-10;
qL = 0;
qH = 10e15;

%At the start eta = 0
eta = 0;
for i = 1:50
    dq = (qL+qH)/2;
    
    if chismall - eta < q/dq
        psi = 1;
    else
        %solve the horrible quadratic
        c = ((a - asmall)/q)*(theta/dtheta)*(-1/chismall)+2*(dq/q)*eta + (dq/q)^2*eta^2 + eta*sigma^2;
        b = ((a-asmall)/q)*(theta/dtheta)*(2/chismall)*(dq/q)*chismall - 2*(dq/q)^2*chismall*eta - chismall*sigma^2;
        a = ((a - asmall)/q)*(theta/dtheta)*(-1/chismall)*(dq/q)^2*(chismall^2);
        %only take positive root
        psi = (-b + sqrt(b^2 - 4*a*c))/(2*a);
    end
    %Try to find sigma_q (You need this to find everything else first
    sigma_q = (sigma)/(1 - (dq/q)*(chismall*psi - eta)) - sigma;
    
    %Try to find sigma_theta as you also need this to find everything else
    sigma_theta = (dtheta/theta)*((chismall*psi - eta)*sigma)/(1 - (dq/q)*(chismall*psi - eta));
    
    %Try to find mu_eta
    mu_eta = (a-iota)/q + ((chi*psi - eta)/(eta))*(sigma + sigma_q)*(-sigma_theta - sigma - sigma_q) + (sigma + sigma_q)*(1-chismall)*(sigma_theta - Sigmasmall);
    
    %Try to find sigma_eta
    sigma_eta = ((chi*psi - eta)/eta)*(sigma + sigma_q);    
    %Try to find mu_q
    mu_q = chismall*(sigma + sigma_q)*(-sigma_theta)- (a-iota)/q + fi(iota) + delta + r;
    %Try to find mu_theta
    mu_theta = rho - r;
    
    
    %Run the ode45
    
end