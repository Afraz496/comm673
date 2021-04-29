clc
clear
close all

color = 'r';

%% Parameters to vary

%Sigma

sigma_vec = [0.25 0.1 0.025];

%Chibar (friction)

chibar = [1 0.5 0.25];

%% Run a baseline
[fout, etaout, dynout] = solve_equilibrium(sigma_vec(1),chibar(1));
N = length(etaout);
% normalize theta, to make sure that theta(eta*) = 1
normalization = fout(N,1);
fout(:,1:2) = fout(:,1:2)/normalization;


figure('Name','Baseline Shooting Method','NumberTitle','off');
subplot(2,2,1); hold on
plot(etaout, fout(:,3), color);   
xlabel('\eta')
ylabel('q');

sigma_eta = dynout(1:N-1,2)./etaout(1:N-1);
subplot(2,2,2); hold on
plot(etaout(1:N-1), dynout(1:N-1,2), color); 
xlabel('\eta')
ylabel('\eta \sigma^{\eta}');

subplot(2,2,3); hold on
plot(etaout(1:N-1), dynout(1:N-1,7), color);  
xlabel('\eta')
ylabel('expert leverage');
axis([0 0.8 0 10]);

subplot(2,2,4); hold on
plot(etaout(1:N-1), dynout(1:N-1,3), color); 
xlabel('\eta')
ylabel('\sigma^q')

%% Prepare the Vectors to show output
dynout_cell = cell(1,3);
fout_cell = cell(1,3);
etaout_cell = cell(1,3);


%% Plots
figure('Name','Sigma','NumberTitle','off')




for i = 1:3
%Call the iterative method on the i-th sigma value
    [fout_cell{i}, etaout_cell{i}, dynout_cell{i}] = solve_equilibrium(sigma_vec(i), chibar(1));
end

subplot(2,2,4); hold on
plot(etaout(1:N-1), dynout(1:N-1,3), color); 


%attempt to plot the i-th sigma value
%Top left
ax1 = nexttile;
plot(ax1, etaout_cell{1}, fout_cell{1}(:,3), 'k');
hold(ax1,'on')
plot(ax1, etaout_cell{2}, fout_cell{2}(:,3),'b');
hold(ax1,'on')
plot(ax1, etaout_cell{3}, fout_cell{3}(:,3),color);
xlim([0 0.8])
xlabel('\eta')
ylabel('q');


sigma_eta1 = dynout_cell{1}(1:N-1,2)./etaout_cell{1}(1:N-1);
N2 = length(etaout_cell{2});
N3 = length(etaout_cell{3});
sigma_eta2 = dynout_cell{2}(1:N2-1,2)./etaout_cell{2}(1:N2-1);
sigma_eta3 = dynout_cell{3}(1:N3-1,2)./etaout_cell{3}(1:N3-1);
ax2 = nexttile;
N = length(etaout_cell{1});
plot(ax2, etaout_cell{1}(1:N-1), sigma_eta1, 'k');
hold(ax2,'on')
N = length(etaout_cell{2});
plot(ax2, etaout_cell{2}(1:N-1), sigma_eta2,'b');
hold(ax2,'on')
N = length(etaout_cell{3});
plot(ax2, etaout_cell{3}(1:N-1), sigma_eta3,color);
hold(ax2,'on')
axis([0 0.8 0 1])
xlabel('\eta')
ylabel('\sigma^{\eta}')


ax3 = nexttile;
N = length(etaout_cell{1});
plot(ax3, etaout_cell{1}(1:N-1), dynout_cell{1}(1:N-1,7), 'k');
hold(ax3, 'on')
N = length(etaout_cell{2});
plot(ax3, etaout_cell{2}(1:N-1), dynout_cell{2}(1:N-1,7),'b');
hold(ax3, 'on')
N = length(etaout_cell{3});
plot(ax3, etaout_cell{3}(1:N-1), dynout_cell{3}(1:N-1,7),color);
hold(ax3, 'on')
axis([0 0.8 0 10]);
xlabel('\eta')
ylabel('expert leverage');

ax4 = nexttile;
N = length(etaout_cell{1});
plot(ax4, etaout_cell{1}(1:N-1), dynout_cell{1}(1:N-1,3), 'k');
hold(ax4,'on')
N = length(etaout_cell{2});
plot(ax4, etaout_cell{2}(1:N-1), dynout_cell{2}(1:N-1,3),'b');
hold(ax4,'on')
N = length(etaout_cell{3});
plot(ax4, etaout_cell{3}(1:N-1), dynout_cell{3}(1:N-1,3),color);
hold(ax4,'on')
xlim([0 0.8])
xlabel('\eta');
ylabel('\sigma^q')


figure('Name','Chibar','NumberTitle','off')


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
xlim([0 0.8])
ylabel('$q$', 'Interpreter', 'latex');
title('\bf{Price of capital good}', 'Interpreter', 'latex')


ax2 = nexttile;
plot(ax2, eta_grid(:,1), sigma_tot_vec(:,1) - sigma_vec(:,1));
hold(ax2,'on')
plot(ax2, eta_grid(:,2), sigma_tot_vec(:,2) - sigma_vec(:,2));
hold(ax2,'on')
plot(ax2, eta_grid(:,3), sigma_tot_vec(:,3) - sigma_vec(:,3));
hold(ax2,'on')
xlim([0 0.8])
ylabel('$\sigma^q$', 'Interpreter', 'latex');
title('\bf{Price volatility}', 'Interpreter', 'latex')


ax3 = nexttile;
plot(ax3, eta_grid(:,1),eta_sigma_eta_vec(:,1));
hold(ax3, 'on')
plot(ax3, eta_grid(:,2),eta_sigma_eta_vec(:,2));
hold(ax3, 'on')
plot(ax3, eta_grid(:,3),eta_sigma_eta_vec(:,3));
hold(ax3, 'on')
xlim([0 0.8])
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