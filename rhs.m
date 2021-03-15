function dx_x = rhs(mu_eta, sigma_eta, mu, eta, x)
%This function converts the second order ODE into a first order ODE using
%the equations given on page 30.
dx_x1 = x(2);
dx_x2 = 2*(mu*x(1) - mu_eta*eta*x(2))/(sigma_eta^2 * eta^2);
dx_x = [dx_x1; dx_x2];
end