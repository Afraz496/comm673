% This file computes the equilibrium in the setting of Brunnermeier and Sannikov
% under the assumptions that (1) experts can issue only debt and not equity and 
% (2) the production set of experts is (a - iota, Phi(iota) - delta) and
% households, (a_ - iota, Phi(iota) - delta_), where the function 
% [Phi iota] = investment(q) solves for Phi and iota
%
% written by Yuliy Sannikov

function [fout, etaout, dynout]= solve_equilibrium(sigmaval,chibarval)
%This is a functional form of the shooting method from Brunnermeri and
%Sannikove (2016). 
%Input: Sigma (volatility), Chibar (skin in the game constraint)
%Outputs: q, sigma_eta, expert_leverage, sigma_q

tic


global qmax

a = 0.11; a_ = 0.05; rho = 0.06; r = 0.05; sigma = sigmaval; delta = 0.03; delta_ = 0.08; chibar = chibarval;
% parameter theta, such that investment iota = Phi + theta Phi^2/2 is
% required to build capital Phi, is defined in the function investment(q)


% SOME PREPARATION: DETERMINE FIRST-BEST Q, AND THE WORST Q (also q(0))
QL = 0; QR = 10; 
for iter = 1:50,
    qmax = (QL + QR)/2;
    [Phi iota] = investment(qmax);
    value = (a - iota)/(r + delta - Phi);
    if iota > a,
        QR = qmax; % price is too large
    elseif value > qmax, % price still too large
        QL = qmax;
    elseif r + delta < Phi,
        'first-best price is infinite' % warn that first-best price is infinite, but allow computation to continue
        QR = qmax;  % without allowing q to grow above the level where Phi > r + delta 
    else
        QR = qmax;
    end
end

% determine q_, the value of capital at eta = 0
QL = 0; QR = 10;
for iter = 1:50,
    q_ = (QL + QR)/2;
    [Phi iota] = investment(q_);
    value = (a_ - iota)/(r + delta_ - Phi);
    if iota > a_,
        QR = q_; % price is too large
    elseif value < q_,  % price still too large
        QR = q_;
    else
        QL = q_;
    end
end


% MAIN PART OF THE CODE:
% solves the system of ODE's starting with boundary conditions F0, i.e.
% theta(0) = 1, theta'(0) = #large negative number# and q(0) = q_; 
% searching for q'(0) such that q'(eta*) and theta'(eta*) reach 0 
% at the same boundary eta*

etaspan = [0 1];
F0 = [1 -1e+10 q_ 0]';  % [theta(0), theta'(0), q(0), q'(0)]

    % note that if theta(eta) satisfies the system of ODE's then so does
    % const*theta(eta), for any positive constant const
    % hence we set theta(0) = 1, and then normalize

options = odeset('RelTol',1e-08,'AbsTol',1e-10, 'events','evntfcn');
        % function 'evntfcn' terminates the integration of the ode 
        % when q exceeds qmax, or f'(eta) or q'(eta) reaches 0 
                                    
odefun = @(eta,f) fnct(eta, f, r, rho, a, a_, delta, delta_, sigma, chibar);

QL = 0; QR = 1e+15;
for iter = 1:50,
    F0(4) = (QL + QR)/2;  % this is q'(0)

    [etaout,fout,TE,YE,IE] = ode45(odefun,etaspan,F0,options);
    
    if IE == 3, % if q'(eta) has reached zero, we 
        QL = F0(4);  % increase q'(0)
    else        % if q(eta) reached qmax or theta'(0) reached 0 
        QR = F0(4);  % we reduce q'(0)
    end

end

% here we are basically done... let me just compute all other variables
% from the same function fnct
N = length(etaout);
dynout = zeros(N, 9);
for n = 1:N,
    [fp dynout(n,:)] = fnct(etaout(n), fout(n,:), r, rho, a, a_, delta, delta_, sigma, chibar);
end





end

