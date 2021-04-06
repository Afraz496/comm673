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
%Find q(0)
%Set Terminal Conditions
f = @(iota)((asmall - iota)/(r - fi(iota) + delta));
g = @(x)-f(x);
f_max = @(iota)((a - iota)/(r - fi(iota) + delta));
g_max = @(x)-f(x);
%this line searches for the maximum value over the specified interval
iota_range = -1/20;
[iota_t,q] = fminsearch(g,iota_range);
q = -1*q;

%The investment function with adjustment costs
fi = @(iota)((log(k*iota + 1))/k);

%ONLY HOLDS IN THE LAST REGION
Sigma = gamma*sigma;

%Solve the equilibrium using equations with 3 parts

%---Part 1---- (The terminal conditions v(n,T) and vsmall(n,T)

%Set the Terminal Conditions for value functions v(eta,T) and vsmall(eta,T)
%Make the grid over eta:
eta = 0.02:0.02:1;
v_t = eta.*(a*eta).^(-gamma); %undefined at eta = 0.
vsmall_t = (1 - eta).*(a*(1-eta)).^(-gamma);
dv_t = (a*eta).^(-gamma) + (-gamma)*a*eta.*((a*eta).^(-gamma-1));
dvsmall_t = -(a*(1-eta).^(-gamma)) + (-gamma)*(-a)*(1-eta).*(-a*(1-eta).^(-gamma -1));
psi = 0.02;
%---Part 2---- (The Static Step)
idx = 1;
qvals = [1];
sigmaqvals = [];
q=1;
delta_eta = 0.02;
for eta = 0.02:delta_eta:0.98
   chi = max([chismall eta]);
   if eta >= chismall
       %This is region 3
       psi = 1;
       sigma_q = 0;
       
       %Find q(eta) from (3.32)
       largepolyroot = 1/k;
       midpolyroot = eta/v_t(idx) + (1-eta)/vsmall_t(idx);
       smallpolyroot = (-(1/k) - a);
       q = roots([largepolyroot midpolyroot smallpolyroot]);
       q = max(q);
   else
       %In region 1
       %Find psi from Equation 3.32
       if psi < 1 || isnan(psi)
            psi = (1/(a-asmall))*(((eta*q)^(1/gamma))/(v_t(idx)^(1/gamma)) + ((1-eta)*q)^(1/gamma)/(vsmall_t(idx))^(1/gamma) - asmall + iota_t);
            
            %Use this psi to find dq
            dqlargepolyroot = ((a-asmall)/q)*((chismall - eta)^2)/q^2;
            dqmidpolyroot = ((a-asmall)/q)*(-2*(chismall*psi - eta)*1/q);
            dqsmallpolyroot = (a-asmall)/q - chismall*(dvsmall_t(idx)/vsmall_t(idx) - dv_t(idx)/v_t(idx) + 1/(eta*(1-eta)))*(chismall*psi - eta)*sigma^2;
            
            dq = roots([dqlargepolyroot dqmidpolyroot dqsmallpolyroot]);
            dq = max(dq);
            q = delta_eta*dq +q;
            
            %using this we can find sigma_q
            sigma_q = sigma/(1 - (chismall - eta)*(dq/q)) - sigma;
            
       end
       
       %Combine equations 3.30, 3.31 and 3.32 to get first order ODE
       
       if psi >= 1
           psi = 1;
       end
       if psi == 1
          % You made it to region 2 
          largepolyroot = 1/k;
          midpolyroot = eta/v_t(idx) + (1-eta)/vsmall_t(idx);
          smallpolyroot = (-(1/k) - a);
          q = roots([largepolyroot midpolyroot smallpolyroot]);
          q = max(q);
          
          %We cheat here to find dq:
          eta_next = eta+delta_eta;
          largepolyroot_next = 1/k;
          midpolyroot_next = eta_next/v_t(idx+1) + (1-eta_next)/vsmall_t(idx+1);
          smallpolyroot_next = (-(1/k) - a);
          q_next = roots([largepolyroot_next midpolyroot_next smallpolyroot_next]);
          q_next = max(q_next);
          
          %Now we can find dq:
          dq = (q_next - q)/delta_eta;
          
          %Now we can find sigma_q
          sigma_q = sigma/(1 - (chismall - eta)*(dq/q)) - sigma;
       end
       
   end
   idx = idx + 1;
   
   %Store values for plotting.
   qvals = [qvals q];
   sigmaqvals = [sigmaqvals sigma_q];
end

figure('Name','Plotting Q','NumberTitle','off')
eta = 0.02:delta_eta:1;
plot(eta, qvals)
xlabel('eta')
ylabel('q')
title('q vs eta')
grid on
%Divide the interval [0,T] into smaller sub intervals as:
delta_t = -0.02;
for t = 1:delta_t:0
    
end


%---Part 3---- (The Time Step)
