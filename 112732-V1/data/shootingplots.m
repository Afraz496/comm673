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

%% Sigma
figure('Name','Sigma','NumberTitle','off')

for i = 1:3
%Call the iterative method on the i-th sigma value
    [fout_cell{i}, etaout_cell{i}, dynout_cell{i}] = solve_equilibrium(sigma_vec(i), chibar(1));
end

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



