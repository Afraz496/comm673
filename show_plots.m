function show_plots(eta_grid, q_vec, sigma_tot_vec, sigma, eta_sigma_eta_vec, eta_mu_eta_vec, N_eta)
% This function generates the standard subplots in the Brunnermeir and
% Sannikov 2016 paper. Plots of q, sigma_q, eta_sigma_eta and eta_mu_eta.
% These figures are generated by the Iterative Method and this function
% allows various plots of different parameter values and acts as a
% placeholder.

subplot(2,2,1);
plot(eta_grid, q_vec);
ylabel('$q$', 'Interpreter', 'latex');
title('\bf{Price of capital good}', 'Interpreter', 'latex')


subplot(2,2,2);
plot(eta_grid, sigma_tot_vec - sigma);
ylabel('$\sigma^q$', 'Interpreter', 'latex');
title('\bf{Price volatility}', 'Interpreter', 'latex')


subplot(2,2,3);
plot(eta_grid,eta_sigma_eta_vec);
ylabel('$\eta\sigma^\eta$', 'Interpreter', 'latex');
title('\bf{Volatility of $\eta$}', 'Interpreter', 'latex')

subplot(2,2,4);
plot(eta_grid,eta_mu_eta_vec, eta_grid, zeros(N_eta,1), 'k:');
ylabel('$\eta\mu^\eta $', 'Interpreter', 'latex');
title('\bf{Drift of $\eta$}', 'Interpreter', 'latex')
end