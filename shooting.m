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
f = @(iota)((asmall - iota)/(r - fi(iota) + delta));
g = @(x)-f(x);
f_max = @(iota)((a - iota)/(r - fi(iota) + delta));
g_max = @(x)-f(x);
%this line searches for the maximum value over the specified interval of
%iota (we need to constrain iota in the range [0,1] since it's a
%reinvestment rate so that can only be 0 -> 1.
iota_range = 0;
[iota,q] = fminsearch(g,iota_range);
theta = 1;
thetainitial = 1;
dtheta = -10e10;
qL = 0;
qH = 10e15;
dthetainitial = -10e10;
%At the start eta = 0
eta = 0.02;
for i = 1:50
    dq = (qL+qH)/2;
    
    if chismall - eta < q/dq
        psi = 1;
    else
        %solve the horrible quadratic
        c = ((a - asmall)/q)*(theta/dtheta)*(-1/chismall)+2*(dq/q)*eta + (dq/q)^2*eta^2 + eta*sigma^2;
        b = ((a-asmall)/q)*(theta/dtheta)*(2/chismall)*(dq/q)*chismall - 2*(dq/q)^2*chismall*eta - chismall*sigma^2;
        aquad = ((a - asmall)/q)*(theta/dtheta)*(-1/chismall)*(dq/q)^2*(chismall^2);
        %only take positive root
        psi = (-b + sqrt(b^2 - 4*aquad*c))/(2*aquad);
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
    [eta_q, q] = ode45(@(eta,q) rhs(mu_eta, sigma_eta, mu_q, eta, q), 0:0.02:0.98, [q dq]);
    [eta_theta, theta] = ode45(@(eta,theta) rhs(mu_eta, sigma_eta, mu_theta, eta, theta), 0:0.02:0.98, [thetainitial dthetainitial]);
    
    dq = q(:,2);
    q = q(:,1);
    dtheta = theta(:,2);
    theta = theta(:,1);
    dq = dq(i);
    q = q(i);
    dtheta = dtheta(i);
    theta = theta(i);
    %If integration terminates for reason (1):
    if q == fminsearch(g_max, 0)
        qH = dq;
    elseif dtheta == 0
        qH = dq;     
    elseif dq == 0
        qL = dq;        
    end
    
    if dtheta == 0 && dq == 0
        eta_optimal = eta_q;
    end
    
    eta = eta + 0.02;
end