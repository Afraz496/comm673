clc
clear
close all

%% Set parameters

% page 40
sigma_vec = [0.1 0.05 0.01];

% page 40
chi_bar_vec = [0.5 0.2 0.1];

% page 41
a_H_vec = [0.03 -0.03 -0.09];

% page 42
gamma_vec = [2 5 0.5];

N_eta = 1000;



%% Invoke the Iterative Method

%check the baseline to see everything is working:
[eta_grid, q_vec, sigma_tot_vec, eta_sigma_eta_vec, eta_mu_eta_vec] = sannikov(sigma_vec(1), a_H_vec(1), chi_bar_vec(1), gamma_vec(1));

%% Plotting figures
figure('Name','Baseline Parameters','NumberTitle','off')
show_plots(eta_grid, q_vec, sigma_tot_vec, sigma_vec(1), eta_sigma_eta_vec, eta_mu_eta_vec, N_eta)

%% Examine Iterative Method for different parameters

%Initilaise returned Sannikov Values in this container
eta_grid = zeros(N_eta, 3);
q_vec = zeros(N_eta, 3);
sigma_tot_vec = zeros(N_eta, 3);
eta_sigma_eta_vec = zeros(N_eta, 3);
eta_mu_eta_vec = zeros(N_eta, 3);

figure('Name','Sigma','NumberTitle','off')




for i = 1:3
%Call the iterative method on the i-th sigma value
    [eta_grid(:,i), q_vec(:,i), sigma_tot_vec(:,i), eta_sigma_eta_vec(:,i), eta_mu_eta_vec(:,i)] = sannikov(sigma_vec(i), a_H_vec(1), chi_bar_vec(1), gamma_vec(1));
end

%attempt to plot the i-th sigma value
%Top left
ax1 = nexttile;
plot(ax1, eta_grid(:,1), q_vec(:,1));
hold(ax1,'on')
plot(ax1, eta_grid(:,2), q_vec(:,2));
hold(ax1,'on')
plot(ax1, eta_grid(:,3), q_vec(:,3));
xlim([0 0.5])
ylabel('$q$', 'Interpreter', 'latex');
title('\bf{Price of capital good}', 'Interpreter', 'latex')


ax2 = nexttile;
plot(ax2, eta_grid(:,1), sigma_tot_vec(:,1) - sigma_vec(:,1));
hold(ax2,'on')
plot(ax2, eta_grid(:,2), sigma_tot_vec(:,2) - sigma_vec(:,2));
hold(ax2,'on')
plot(ax2, eta_grid(:,3), sigma_tot_vec(:,3) - sigma_vec(:,3));
hold(ax2,'on')
xlim([0 0.5])
ylabel('$\sigma^q$', 'Interpreter', 'latex');
title('\bf{Price volatility}', 'Interpreter', 'latex')


ax3 = nexttile;
plot(ax3, eta_grid(:,1),eta_sigma_eta_vec(:,1));
hold(ax3, 'on')
plot(ax3, eta_grid(:,2),eta_sigma_eta_vec(:,2));
hold(ax3, 'on')
plot(ax3, eta_grid(:,3),eta_sigma_eta_vec(:,3));
hold(ax3, 'on')
xlim([0 0.5])
ylabel('$\eta\sigma^\eta$', 'Interpreter', 'latex');
title('\bf{Volatility of $\eta$}', 'Interpreter', 'latex')

ax4 = nexttile;
plot(ax4, eta_grid(:,1),eta_mu_eta_vec(:,1), eta_grid(:,1), zeros(N_eta,1), 'k:');
hold(ax4,'on')
plot(ax4, eta_grid(:,2),eta_mu_eta_vec(:,2),'r', eta_grid(:,2), zeros(N_eta,1), 'k:');
hold(ax4,'on')
plot(ax4, eta_grid(:,3),eta_mu_eta_vec(:,3),'y', eta_grid(:,3), zeros(N_eta,1), 'k:');
hold(ax4,'on')
xlim([0 0.5])
ylabel('$\eta\mu^\eta $', 'Interpreter', 'latex');
title('\bf{Drift of $\eta$}', 'Interpreter', 'latex')


figure('Name','chi_bar','NumberTitle','off')

for i = 1:3
%Call the iterative method on the i-th chibar value
    [eta_grid(:,i), q_vec(:,i), sigma_tot_vec(:,i), eta_sigma_eta_vec(:,i), eta_mu_eta_vec(:,i)] = sannikov(sigma_vec(1), a_H_vec(1), chi_bar_vec(i), gamma_vec(1));
end

%attempt to plot the i-th chibar value
%Top left
ax1 = nexttile;
plot(ax1, eta_grid(:,1), q_vec(:,1));
hold(ax1,'on')
plot(ax1, eta_grid(:,2), q_vec(:,2));
hold(ax1,'on')
plot(ax1, eta_grid(:,3), q_vec(:,3));
xlim([0 0.5])
ylabel('$q$', 'Interpreter', 'latex');
title('\bf{Price of capital good}', 'Interpreter', 'latex')


ax2 = nexttile;
plot(ax2, eta_grid(:,1), sigma_tot_vec(:,1) - sigma_vec(:,1));
hold(ax2,'on')
plot(ax2, eta_grid(:,2), sigma_tot_vec(:,2) - sigma_vec(:,1));
hold(ax2,'on')
plot(ax2, eta_grid(:,3), sigma_tot_vec(:,3) - sigma_vec(:,1));
hold(ax2,'on')
xlim([0 0.5])
ylabel('$\sigma^q$', 'Interpreter', 'latex');
title('\bf{Price volatility}', 'Interpreter', 'latex')


ax3 = nexttile;
plot(ax3, eta_grid(:,1),eta_sigma_eta_vec(:,1));
hold(ax3, 'on')
plot(ax3, eta_grid(:,2),eta_sigma_eta_vec(:,2));
hold(ax3, 'on')
plot(ax3, eta_grid(:,3),eta_sigma_eta_vec(:,3));
hold(ax3, 'on')
xlim([0 0.5])
ylabel('$\eta\sigma^\eta$', 'Interpreter', 'latex');
title('\bf{Volatility of $\eta$}', 'Interpreter', 'latex')

ax4 = nexttile;
plot(ax4, eta_grid(:,1),eta_mu_eta_vec(:,1), eta_grid(:,1), zeros(N_eta,1), 'k:');
hold(ax4,'on')
plot(ax4, eta_grid(:,2),eta_mu_eta_vec(:,2),'r', eta_grid(:,2), zeros(N_eta,1), 'k:');
hold(ax4,'on')
plot(ax4, eta_grid(:,3),eta_mu_eta_vec(:,3),'y', eta_grid(:,3), zeros(N_eta,1), 'k:');
hold(ax4,'on')
xlim([0 0.5])
ylabel('$\eta\mu^\eta $', 'Interpreter', 'latex');
title('\bf{Drift of $\eta$}', 'Interpreter', 'latex')

figure('Name','a_H','NumberTitle','off')
for i = 1:3
%Call the iterative method on the i-th sigma value
    i
    [eta_grid(:,i), q_vec(:,i), sigma_tot_vec(:,i), eta_sigma_eta_vec(:,i), eta_mu_eta_vec(:,i)] = sannikov(sigma_vec(1), a_H_vec(i), chi_bar_vec(1), gamma_vec(1));
end

%attempt to plot the i-th a_H value
%Top left
ax1 = nexttile;
plot(ax1, eta_grid(:,1), q_vec(:,1));
hold(ax1,'on')
plot(ax1, eta_grid(:,2), q_vec(:,2));
hold(ax1,'on')
plot(ax1, eta_grid(:,3), q_vec(:,3));
ylabel('$q$', 'Interpreter', 'latex');
title('\bf{Price of capital good}', 'Interpreter', 'latex')


ax2 = nexttile;
plot(ax2, eta_grid(:,1), sigma_tot_vec(:,1) - sigma_vec(:,1));
hold(ax2,'on')
plot(ax2, eta_grid(:,2), sigma_tot_vec(:,2) - sigma_vec(:,1));
hold(ax2,'on')
plot(ax2, eta_grid(:,3), sigma_tot_vec(:,3) - sigma_vec(:,1));
hold(ax2,'on')
ylabel('$\sigma^q$', 'Interpreter', 'latex');
title('\bf{Price volatility}', 'Interpreter', 'latex')


ax3 = nexttile;
plot(ax3, eta_grid(:,1),eta_sigma_eta_vec(:,1));
hold(ax3, 'on')
plot(ax3, eta_grid(:,2),eta_sigma_eta_vec(:,2));
hold(ax3, 'on')
plot(ax3, eta_grid(:,3),eta_sigma_eta_vec(:,3));
hold(ax3, 'on')
ylabel('$\eta\sigma^\eta$', 'Interpreter', 'latex');
title('\bf{Volatility of $\eta$}', 'Interpreter', 'latex')

ax4 = nexttile;
plot(ax4, eta_grid(:,1),eta_mu_eta_vec(:,1), eta_grid(:,1), zeros(N_eta,1), 'k:');
hold(ax4,'on')
plot(ax4, eta_grid(:,2),eta_mu_eta_vec(:,2),'r', eta_grid(:,2), zeros(N_eta,1), 'k:');
hold(ax4,'on')
plot(ax4, eta_grid(:,3),eta_mu_eta_vec(:,3),'y', eta_grid(:,3), zeros(N_eta,1), 'k:');
hold(ax4,'on')
ylabel('$\eta\mu^\eta $', 'Interpreter', 'latex');
title('\bf{Drift of $\eta$}', 'Interpreter', 'latex')


figure('Name','gamma','NumberTitle','off')
for i = 1:3
%Call the iterative method on the i-th sigma value
    [eta_grid(:,i), q_vec(:,i), sigma_tot_vec(:,i), eta_sigma_eta_vec(:,i), eta_mu_eta_vec(:,i)] = sannikov(sigma_vec(1), a_H_vec(1), chi_bar_vec(1), gamma_vec(i));
    show_plots(eta_grid(:,i), q_vec(:,i), sigma_tot_vec(:,i), sigma_vec(1), eta_sigma_eta_vec(:,i), eta_mu_eta_vec(:,i), N_eta)
    hold on
end

%attempt to plot the i-th gamma value
%Top left
ax1 = nexttile;
plot(ax1, eta_grid(:,1), q_vec(:,1));
hold(ax1,'on')
plot(ax1, eta_grid(:,2), q_vec(:,2));
hold(ax1,'on')
plot(ax1, eta_grid(:,3), q_vec(:,3));
ylabel('$q$', 'Interpreter', 'latex');
title('\bf{Price of capital good}', 'Interpreter', 'latex')


ax2 = nexttile;
plot(ax2, eta_grid(:,1), sigma_tot_vec(:,1) - sigma_vec(:,1));
hold(ax2,'on')
plot(ax2, eta_grid(:,2), sigma_tot_vec(:,2) - sigma_vec(:,1));
hold(ax2,'on')
plot(ax2, eta_grid(:,3), sigma_tot_vec(:,3) - sigma_vec(:,1));
hold(ax2,'on')
ylabel('$\sigma^q$', 'Interpreter', 'latex');
title('\bf{Price volatility}', 'Interpreter', 'latex')


ax3 = nexttile;
plot(ax3, eta_grid(:,1),eta_sigma_eta_vec(:,1));
hold(ax3, 'on')
plot(ax3, eta_grid(:,2),eta_sigma_eta_vec(:,2));
hold(ax3, 'on')
plot(ax3, eta_grid(:,3),eta_sigma_eta_vec(:,3));
hold(ax3, 'on')
ylabel('$\eta\sigma^\eta$', 'Interpreter', 'latex');
title('\bf{Volatility of $\eta$}', 'Interpreter', 'latex')

ax4 = nexttile;
plot(ax4, eta_grid(:,1),eta_mu_eta_vec(:,1), eta_grid(:,1), zeros(N_eta,1), 'k:');
hold(ax4,'on')
plot(ax4, eta_grid(:,2),eta_mu_eta_vec(:,2), eta_grid(:,2), zeros(N_eta,1), 'k:');
hold(ax4,'on')
plot(ax4, eta_grid(:,3),eta_mu_eta_vec(:,3), eta_grid(:,3), zeros(N_eta,1), 'k:');
hold(ax4,'on')
ylabel('$\eta\mu^\eta $', 'Interpreter', 'latex');
title('\bf{Drift of $\eta$}', 'Interpreter', 'latex')