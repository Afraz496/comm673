clc
clear
close all

%Baseline Parameters
rho = 0.06;
r = 0.05;
a = 0.11;
asmall = 0.03;
delta = 0.05;
sigma = 0.1;
chismall = 0.5;
gamma = 2;
k = 10;

%Solve the equilibrium using equations with 3 parts

%---Part 1---- (The terminal conditions v(n,T) and vsmall(n,T)

%Set the Terminal Conditions for value functions v(eta,T) and vsmall(eta,T)
%Make the grid over eta:
eta = 0:0.01:1;
v_t = eta*(a*eta)^(-1*gamma);
vsmall_t = (1 - eta)*(a(1-eta))^(-1*gamma);

%Divide the interval [0,T] into smaller sub intervals as:
delta_t = 0.01;
t = 0:delta_t:1;

%---Part 2---- (The Static Step)

%---Part 3---- (The Time Step)
