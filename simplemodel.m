clc
clear
close all

%Simple Real Model Plots

%Parameter Values
a = 0.11;
rho = 0.05;
sigma = 0.1;
K = 10;
delta = 0.0336; %Reverse Engineered this from the graph on the paper (it is exactly fi(iota)
fi = @(iota)((log(K*iota + 1))/K);


%Generate Plots for Figure 2 in the paper (page 16)

% Top Left Panel - q vs eta (q: price of capital, eta: Share of experts
% wealth in [0,1])

figure('Name', 'Top-Left','NumberTitle','off')

eta = 0:0.02:1;
q = 1.4*ones(length(eta),1);
plot(eta, q)
axis([0 1 0 1.8])
xlabel('$\eta$','Interpreter','latex')
ylabel('q')

% Top Right Panel - risk free rate vs eta

% This is equation 2.9 from page 14
iota = (q - 1)./K;
rf = rho*ones(length(iota),1) + fi(iota) - delta*ones(length(iota),1) - (sigma^2./eta)';

figure('Name','Top-Right','NumberTitle','off')
plot(eta, rf)
axis([0 1 -inf 0.1])
xlabel('$\eta$','Interpreter','latex')
ylabel('$r_t$, risk-free rate','Interpreter','latex')


%Bottom Left Panel
etasigma = normpdf(eta,0.5,1);
figure('Name','Bottom-Left','NumberTitle','off')
plot(eta,eta.*etasigma)
