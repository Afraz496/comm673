% ***********************************************************************
% This file computes the equilibrium for the model in Section 3 of the
% Brunnermeier and Sannikov's handbook chapter:
% "Money, Macro, and Finance, A Continuous-Time Approach" (2016).
%
% Assumptions:
% (1) Experts can issue debt and equity, but must retain at least fraction alpha*psi of risk and
% (2) Both experts and households have CRRA utility with relative risk aversion gamma
% The investment function is assumed to be of the form Phi(iota) =
% log(kappa*iota + 1)/kappa
%
% Adapted from the original code of Yuliy Sannikov
% ***********************************************************************

function [eta_grid, q_vec, sigma_tot_vec, eta_sigma_eta_vec, eta_mu_eta_vec] = sannikov(sigma, a_H, chi_bar, gamma)

Newton_sw = 1;  % Do you want to use Newton method for nonlinear euqation solution? 1=yes, 0=no (use fsolve)
stationary_distr_sw = 0;   %=0: don't compute the stationary distribution; =1: compute the stationary distribution;

tic
%% Parameters

%Preferences
gamma_E = gamma;
gamma_H = gamma;

% Technology
kappa = 10;
delta = 0.03;
rho_E = 0.06;
rho_H = 0.05;
a_E   = 0.11;



%% Grid for net worth share of expert: eta
% uneven grid on [0, 1] with more points near 0 and 1

N_T      = 40;
N_eta    = 1000;
eta_low  = 2.998000000000000e-06;
eta_high = 0.999997002000000;

eta_transf = @(x)  3*x.^2 - 2*x.^3;

z_low  = fzero(@(x) eta_transf(x) - eta_low,  1e-4);
z_high = fzero(@(x) eta_transf(x) - eta_high, 1-1e-4);

if z_low<0 || z_high>1 || z_low>1 || z_high<0
    fprintf('z_low = %f\n' ,z_low)
    fprintf('z_high = %f\n',z_high)
    error('Eta map not correct!')
end



z_grid = linspace(z_low, z_high, N_eta)';

eta_grid = eta_transf(z_grid);

if eta_grid(1)<0 || eta_grid(end)>1
    error('eta<0 or eta>1!')
end

d_eta = diff(eta_grid);



%% Terminal conditions for value functions v_E and v_H  (C_T = a*K*eta: consume everything at time T)

v_E  = a_E^(-gamma_E) * eta_grid.^(1 - gamma_E);
v_H  = a_H^(-gamma_H) * (1-eta_grid).^(1 - gamma_H);  %Changed to a_H and works fine (actually fixes the code)

q_vec           = NaN(N_eta,1);
q_prime_vec     = NaN(N_eta,1);
sigma_tot_vec   = NaN(N_eta,1);
psi_vec         = NaN(N_eta,1);
eta_mu_eta_vec  = zeros(N_eta,1);
eta_sig_eta_vec = zeros(N_eta,1);
mu_v_E_vec      = NaN(N_eta,1);
mu_v_H_vec      = NaN(N_eta,1);
sigma_v_E_vec   = NaN(N_eta,1);
sigma_v_H_vec   = NaN(N_eta,1);

payoff_vec      = zeros(N_eta,1); %Used in the PDE function


psi_vec(1)       = 0;
sigma_tot_vec(1)   = 0;
sigma_v_E_vec(1)   = 0;
sigma_v_H_vec(1)   = 0;


chi_vec =  max(chi_bar,eta_grid);


for t = 1:N_T
    disp(t)
    %% Initialization of term not depending on q, Psi, Chi, and sigma_tot
    
    CK_no_q_E = (eta_grid./v_E).^(1/gamma_E);
    CK_no_q_H = ((1-eta_grid)./v_H).^(1/gamma_H);
    CK_no_q   = CK_no_q_E + CK_no_q_H;
    
    
    v_prime_E  = diff(v_E)./d_eta;
    v_prime_H  = diff(v_H)./d_eta;
    
    v_prime_over_v_E =  v_prime_E./v_E(2:end);
    v_prime_over_v_H =  v_prime_H./v_H(2:end);
    
    
    FOC_psi_term = -v_prime_over_v_E + v_prime_over_v_H +...
        1./((1-eta_grid(2:end)).*eta_grid(2:end));
    
    
    %% Initial condition [q,psi,sigma_tot] at eta=0 (no experts)
    
    
    % Find the market clearing price of capital q(0)
    qL  = 0;
    qR  = a_H*kappa + 1;
    
    discr = 1;
    
    if a_H == -0.09 || gamma == 0.5
        discr_condition = 1e-100;
    else
        discr_condition = 1e-10;
    end
    
    %original discr condition was 1e-10
    while abs(discr) > discr_condition
        %     for j = 1:30                              %NOTE: Change this after confirming solution
        q          = (qL + qR)/2;
        iota       = (q - 1)/kappa;
        CK_0       = a_H - iota;
        
        log(CK_0)
        discr = log(q)/gamma_H + log(CK_no_q(1)) - log(CK_0); % mkt clearing
        if discr>0
            qR = q;
        else
            qL = q;
        end
        %         fprintf('q = %f, discr = %f\n', q, discr);
    end
    
    q_vec(1) = q;
    q_old    = q;
    
    psi       = psi_vec(1);
    sigma_tot = sigma;
    
    
    %% Newton method
    %find q, psi and sigma_tot = sigma + sigma^q
    
    for i_eta = 2:N_eta
        
        eta = eta_grid(i_eta);
        
        % ---- STATIC STEP ------
        if Newton_sw == 1  %Use Newton method to solve the nonlinear system of equation
            
            
            x0 = [q, psi, sigma_tot];
            ER = static_equil(x0, q_old, eta, d_eta(i_eta-1), CK_no_q_E(i_eta), CK_no_q_H(i_eta),...
                chi_bar, sigma, gamma_E, gamma_H, a_E, a_H, kappa, FOC_psi_term(i_eta-1));
            
            dER = compute_jacobian(q, q_old, psi, sigma_tot, eta, d_eta(i_eta-1),...
                gamma_E,gamma_H, CK_no_q_E(i_eta), CK_no_q_H(i_eta),  kappa, a_E, a_H,...
                sigma, chi_bar, FOC_psi_term(i_eta-1));
            
            new_guess = [q; psi; sigma_tot] - dER\ER;
            
            if sum(new_guess<0)>0
                t
                i_eta
                new_guess
                error('Solution to static step is wrong!')
            end
            
        elseif Newton_sw==0  %Use non-linear equation solver
            x0 = [q, psi, sigma_tot];
            [x, fval, exitflag] = fsolve(@(x) static_equil(x, q_old, eta, d_eta(i_eta-1), CK_no_q_E(i_eta), CK_no_q_H(i_eta),...
                chi_bar, sigma, gamma_E, gamma_H, a_E, a_H, kappa, FOC_psi_term(i_eta-1)),x0,...
                optimoptions('fsolve','Display','off'));
            %  x
            %  fval'
            %  exitflag
            new_guess = x;
            
        end
        
        % If the boundary of the crisis regime has been reached, we have
        % psi = 1 from now on
        
        if new_guess(2) > 1  %psi>1
            break
        end
        
        % update new guesses at current eta
        q         = new_guess(1);
        psi       = new_guess(2);
        sigma_tot = new_guess(3);
        
        
        % Collect results
        
        q_vec(i_eta)   = new_guess(1);
        psi_vec(i_eta) = new_guess(2);
        sigma_tot_vec(i_eta) = new_guess(3);
        q_prime_vec(i_eta)  = (q_vec(i_eta) - q_old)/d_eta(i_eta-1);
        
        q_old  = q;  % For use in derivatives of q w.r.t. eta
        
    end
    
    %% Newton method for psi = 1, for remaining eta
    
    i_eta_psi_1 = i_eta;
    
    for i_eta = i_eta_psi_1:N_eta
        psi = 1;
        eta = eta_grid(i_eta); 
        
        
        
        if Newton_sw == 1  %Use Newton method to solve the nonlinear equation
            
            ER = static_equil_psi_1(q, psi, gamma_E, gamma_H,CK_no_q_E(i_eta), CK_no_q_H(i_eta), kappa, a_E, a_H);
            
            iota = (q - 1)/kappa;
            dER_dq = 1/(q*gamma_E)  + 1/(a_E - iota) * 1/kappa;
            new_guess = q - ER/dER_dq;
            
        else %Use non-linear equation solver
            x0 = q;
            [x, fval] = fzero(@(x) static_equil_psi_1(x, psi, gamma_E, gamma_H,CK_no_q_E(i_eta), CK_no_q_H(i_eta), kappa, a_E, a_H), x0);
            new_guess = x;
        end
        
        q      =  new_guess;
        q_vec(i_eta)   = q;
        psi_vec(i_eta) = psi;
        
        q_prime_vec(i_eta)  = (q_vec(i_eta) - q_vec(i_eta-1))/d_eta(i_eta-1);
        sigma_tot_vec(i_eta) = 1./(1 - (chi_vec(i_eta) - eta)*q_prime_vec(i_eta)/q_vec(i_eta))*sigma;
        
    end
    
    
    
    %% Computing the PDE
    %Volatility of Eta
    
    eta_sigma_eta_vec    = (psi_vec.*chi_vec - eta_grid).*sigma_tot_vec;
    %   sigma_eta_vec(N_eta) = 0;  %Redundant
    iota_vec        = (q_vec - 1)/kappa;
    Phi_vec         = log(q_vec)/kappa;   %Recall: adjustment cost function: phi(i) = 1/kappa*log(kappa*iota + 1)
    A_psi_vec       = a_E*psi_vec + a_H*(1 - psi_vec);
    CK_vec          = A_psi_vec - iota_vec;
    
    
    % for the law of motion of Eta, C/N = (Eta Q)^(1/gamma - 1)/V^(1/gamma)
    % recall that GS1 = (Eta./V).^(1/gamma);
    
    CN_E_vec  = CK_no_q_E.*q_vec.^(1/gamma_E - 1)./eta_grid;
    CN_H_vec  = CK_no_q_H.*q_vec.^(1/gamma_H - 1)./(1-eta_grid);
    
    sigma_v_E_vec(2:N_eta) = v_prime_over_v_E.*eta_sigma_eta_vec(2:N_eta);
    sigma_v_H_vec(2:N_eta) = v_prime_over_v_H.*eta_sigma_eta_vec(2:N_eta);
    
    riskPrice_E_vec = -sigma_v_E_vec + eta_sigma_eta_vec./eta_grid     + sigma_tot_vec + (gamma_E-1)*sigma;
    %     riskPrice_H_vec = -sigma_v_H_vec + eta_sigma_eta_vec./(1-eta_grid) + sigma_tot_vec + (gamma_H-1)*sigma;  %THIS IS WRONG!!!!!
    riskPrice_H_vec = -sigma_v_H_vec - eta_sigma_eta_vec./(1-eta_grid) + sigma_tot_vec + (gamma_H-1)*sigma;  %THIS IS THE CORRECT VERSION!!!!!
    
    
    eta_mu_eta_vec(2:N_eta-1) = eta_grid(2:N_eta-1) .* ( (a_E - iota_vec(2:N_eta-1))./q_vec(2:N_eta-1) - CN_E_vec(2:N_eta-1) ) + ...
        eta_sigma_eta_vec(2:N_eta-1) .* (riskPrice_E_vec(2:N_eta-1) - sigma_tot_vec(2:N_eta-1)) + ...
        eta_grid(2:N_eta-1) * (1-chi_bar) .* (riskPrice_E_vec(2:N_eta-1) - riskPrice_H_vec(2:N_eta-1)) .* sigma_tot_vec(2:N_eta-1); %NOTE: need to change to adapt to the slides version of the model!!!!!!!!!!!!!!!
    
    
    mu_v_E_vec = rho_E - CN_E_vec - (1 - gamma_E)*(Phi_vec - delta - 1/2*gamma_E*sigma^2 + sigma*sigma_v_E_vec);
    mu_v_H_vec = rho_H - CN_H_vec - (1 - gamma_H)*(Phi_vec - delta - 1/2*gamma_H*sigma^2 + sigma*sigma_v_H_vec);
    
    
    
    %% Updating V and V_
    % The parameter lambda = dt/(1 + dt). lambda = 1 corresponds to policy iteration
    % it is more agressive to set lambda closer to 1, but the code may not converge
    lambda = 0.8;
    
    v_E  = solve_PDE_implicit(eta_grid, mu_v_E_vec, eta_mu_eta_vec, eta_sigma_eta_vec, payoff_vec, v_E, lambda);
    v_H  = solve_PDE_implicit(eta_grid, mu_v_H_vec, eta_mu_eta_vec, eta_sigma_eta_vec, payoff_vec, v_H, lambda);
    
end

toc

%% Determining mu_q and the riskfree rate [Not sure about this]


sigma_N_E_vec = (chi_vec.*psi_vec)./eta_grid.*sigma_tot_vec;
sigma_N_H_vec = (1-chi_vec.*psi_vec)./(1-eta_grid).*sigma_tot_vec;


q_second_vec      = NaN(N_eta,1);
q_second_vec(3:N_eta) = diff(q_prime_vec(2:end))./d_eta(2:end);

% mu_q_vec = NaN(N_eta,1);
q_mu_q_vec = q_prime_vec.*eta_mu_eta_vec +...
    1/2*q_second_vec.*(eta_sigma_eta_vec).^2;


mu_eta_vec    = eta_mu_eta_vec./eta_grid;
sigma_eta_vec = eta_sigma_eta_vec./eta_grid;
mu_q_vec      = q_mu_q_vec./q_vec;
sigma_q_vec   = sigma_tot_vec - sigma;
mu_qK_vec     = mu_q_vec + Phi_vec - delta + sigma_q_vec*sigma;
sigma_qK_vec  = sigma_tot_vec;
mu_N_E_vec    = mu_eta_vec + mu_qK_vec + sigma_eta_vec.*sigma_qK_vec;
rf_vec_v1    = mu_N_E_vec - (chi_vec.*psi_vec)./eta_grid.*sigma_tot_vec;


rf_vec_E = (a_E - iota_vec)./q_vec + Phi_vec - delta + mu_q_vec + sigma_tot_vec -...
    (riskPrice_E_vec.*chi_vec + riskPrice_H_vec.*(1-chi_vec)).*sigma_tot_vec;

rf_vec_H = (a_H - iota_vec)./q_vec + Phi_vec - delta + mu_q_vec + sigma_tot_vec -...
    riskPrice_H_vec.*sigma_tot_vec;


%% Checking Equilibrium Conditions

FOC_1 = log(CK_vec) - log(A_psi_vec - iota_vec);
FOC_2 = sigma_tot_vec(2:end).*(q_vec(2:end) - (q_vec(2:end) - q_vec(1:end-1)).*(chi_bar*psi_vec(2:end) - eta_grid(2:end))./d_eta) - sigma*q_vec(2:end);
FOC_3 = a_E - a_H -  q_vec(2:end).*chi_bar.*(chi_bar*psi_vec(2:end) - eta_grid(2:end)).*sigma_tot_vec(2:end).^2.*FOC_psi_term;




%% Stationary distribution (Kolmogorov Forward Equation)

if stationary_distr_sw == 1
    [f, A] = compute_stationary_distribution(eta_grid, eta_mu_eta_vec, eta_sigma_eta_vec);
end





%% AUXILIARY FUNCTIONS

function [CK, iota, A_psi] = CK_and_investment(q, psi, CK_no_q_E, CK_no_q_H, gamma_E, gamma_H, a_E, a_H, kappa)

CK_E   = q^(1/gamma_E) * CK_no_q_E;
CK_H   = q^(1/gamma_H) * CK_no_q_H;
CK     = CK_E + CK_H;

iota  = (q - 1)/kappa;
A_psi = a_E*psi + a_H*(1 - psi);


end


function dER = compute_jacobian(q, q_old, psi, sigma_tot, eta, d_eta, ...
    gamma_E,gamma_H, CK_no_q_E, CK_no_q_H, kappa, a_E, a_H, ...
    sigma, chi_bar, FOC_psi_term)


iota  = (q - 1)/kappa;
A_psi = a_E*psi + a_H*(1 - psi);

q_prime   = (q - q_old)/d_eta;


% matrix of derivatives of errors
% Note that the jacobian is computed at q(i_eta) = q_old = q(i_eta-1), hence q_prime(i_eta) = 0!
% could shorten it since q_old = q


dER = NaN(3,3);

dER1_dq = 1/(q^(1/gamma_E)*CK_no_q_E + q^(1/gamma_H)*CK_no_q_H)*...
    (1/gamma_E*q^(1/gamma_E-1)*CK_no_q_E + ...
    1/gamma_H*q^(1/gamma_H-1)*CK_no_q_H)+ ...
    1/(A_psi - iota)*1/kappa;
dER1_dpsi = -(a_E - a_H)/(A_psi-iota);
dER1_dsigma_tot = 0;


dER2_dq   = sigma_tot*(1 - (chi_bar*psi - eta)/d_eta) - sigma;
dER2_dpsi = -sigma_tot * q_prime * chi_bar;
dER2_dsigma_tot = q - q_prime*(chi_bar*psi - eta);


dER3_dq   = - chi_bar*FOC_psi_term*(chi_bar*psi - eta)*sigma_tot^2;
dER3_dpsi = -q*chi_bar^2*sigma_tot^2*FOC_psi_term;
dER3_dsigma_tot = -2*q*chi_bar*(chi_bar*psi - eta)*sigma_tot*FOC_psi_term;



dER(1,:) = [dER1_dq, dER1_dpsi, dER1_dsigma_tot];
dER(2,:) = [dER2_dq, dER2_dpsi, dER2_dsigma_tot];
dER(3,:) = [dER3_dq, dER3_dpsi, dER3_dsigma_tot];


end


function  static_equil_val = static_equil(x, q_old, eta, d_eta, CK_no_q_E, CK_no_q_H,...
    chi_bar, sigma, gamma_E, gamma_H, a_E, a_H, kappa, FOC_psi_term)
q         = x(1);
psi       = x(2);
sigma_tot = x(3);

[CK, iota, A_psi] = CK_and_investment(q, psi, CK_no_q_E, CK_no_q_H, gamma_E, gamma_H, a_E, a_H, kappa);

q_prime = (q-q_old)/d_eta;

% Equilibrium conditions
Mkt_clearing  = log(CK) - log(A_psi - iota);
sigma_q_equil = sigma_tot*(q - q_prime*(chi_bar*psi - eta)) - sigma*q;
FOC_psi       =  (a_E - a_H) - q * chi_bar*FOC_psi_term*(chi_bar*psi - eta)*sigma_tot^2;

static_equil_val = [Mkt_clearing; sigma_q_equil; FOC_psi];


end


function  static_equil_psi_1_val = static_equil_psi_1(x, psi, gamma_E, gamma_H,...
    CK_no_q_E, CK_no_q_H, kappa, a_E, a_H)
q  = x;

[CK, iota, A_psi] = CK_and_investment(q, psi, CK_no_q_E, CK_no_q_H, gamma_E, gamma_H, a_E, a_H, kappa);

% Equilibrium conditions
Mkt_clearing = log(CK) - log(A_psi - iota);

static_equil_psi_1_val = Mkt_clearing;

end


function F = solve_PDE_implicit(X, R, MU, S, G, V, lambda)
% F = solve_PDE_implicit(X, R, MU, S, G, V, lambda)
% where X, R, MU, S, and G are column vectors of the same length
%
% X = [X(1), X(2) ... X(N)]' is the state space (an increasing grid)
% R is the discount rate minus growth
% MU is a drift vector of length N, with MU(1) >= 0 and MU(N) <= 0
% S is a volatility vector of lenght N, with S(1) = S(N) = 0
% G is a payoff flow of length N
%
% we are solving the value function equation
% R(X)*F(t,X) = G(X) + MU(X)*F_x(t,X) + S(X)^2*F_xx(t,X)/2 + F_t(t,X)
% where F_t(t,X) = (F(t,X) - F(t-dt))/dt, with F(t,X) = V(t,X) known
%
% to solve the equation over a small time interval dt with continuation
% value function V, set lambda = dt/(1 + dt) (or simply lambda = dt)
%
% to get the stationary solution of this equation, set lambda = 1:
% then V has no effect on the solution
%
% Adapted from the original code of Yuliy Sannikov


if or(MU(1) < 0, MU(end) > 0)
    MU(1)
    MU(end)
    error('Error: not true that MU(1) >= 0 and MU(N) <= 0');
end

if or(S(1) ~= 0, S(end) ~= 0)
    S(1)
    S(end)
    error('Error: not true that S(1) and S(N) both nonzero');
end

N = length(X);
dX = X(2:N) - X(1:N-1);


%-------- Construction of the sparse tri-diagonal matrix -------

S2_dX2 = [0;S(2:N-1).^2 ./ (dX(1:N-2) + dX(2:N-1));0];

uD = lambda*(-max(MU(1:N-1),0) - S2_dX2(1:N-1))./dX;
lD = lambda*(min(MU(2:N),0) - S2_dX2(2:N))./dX;
mD = lambda*R + 1-lambda - [uD;0] - [0;lD];


upper_diag = [0;uD];
lower_diag = [lD;0];
main_diag  = mD;

A = spdiags([lower_diag, main_diag, upper_diag], -1:1, N,N );

F = A\(G*lambda + V*(1 - lambda));



%------------------- Sannikov's version ------------------------------
% S0 = zeros(N,1); S0(2:N-1) = S(2:N-1).^2./(dX(1:N-2) + dX(2:N-1));
% DU = zeros(N,1); DU(2:N) = - (max(MU(1:N-1),0) + S0(1:N-1))./dX*lambda;
% DD = zeros(N-1,1); DD = - (max(-MU(2:N),0) + S0(2:N))./dX*lambda;
%
% D0 = (1 - lambda)*ones(N,1) + lambda*R;
% D0(1:N-1) = D0(1:N-1) - DU(2:N); D0(2:N) = D0(2:N) - DD;
% A = spdiags(D0,0,N,N) + spdiags(DU,1,N,N) + spdiags(DD(1:N-1),-1,N,N);
% F = A\(G*lambda + V*(1 - lambda));



end

function [f, A] = compute_stationary_distribution(X, MU, S)

N   = length(X);
dX  = X(2:N) - X(1:N-1);
dX2 = dX(1:N-2) + dX(2:N-1);
S2  = S.^2;

%-------- Construction of the sparse tri-diagonal matrix -------



uD =  (-min(MU(2:N),0)  + [S2(2:N-1)./dX2; 0])./dX;
lD =  (max(MU(1:N-1),0) + [0; S2(2:N-1)./dX2])./dX;
mD =  -[0;uD] - [lD;0];


upper_diag = [0;uD];
lower_diag = [lD;0];
main_diag  = mD;

A = spdiags([lower_diag, main_diag, upper_diag], -1:1, N,N );

%compute largest real eigenvalues and associated eigenvector of matrix A
% [eigenVec,eigenvalues] = eigs(A',1,'largestreal');  %From https://benjaminmoll.com/wp-content/uploads/2020/02/HACT_Numerical_Appendix.pdf footnote 5, p.9

% f = eigenVec;  %%%% NOT SURE ABOUT THIS!!!!



% f_old = [0;ones(N-2,1)*1/(N-2);0];

n_zeros = 10;
f_old = [zeros(n_zeros,1);ones(N-2*n_zeros,1)*1/(N-2*n_zeros);zeros(n_zeros,1)];



dt = min(abs(dX))/max(abs(MU))/10;
% dt = 0.0001;

discr = 1;
iter = 0;
while discr>1e-9
    % for i=1:10000
    
    iter = iter + 1;
    f_new = (A*dt + speye(N))*f_old;
    
    discr = sum(abs(f_new-f_old));
    if mod(iter,10000)==0
        fprintf('iter %d, discr = %g\n', iter, discr);
    end
    f_old = f_new;
end

f = f_new/(f_new(2:end)'*dX);  %normalization

end
end

