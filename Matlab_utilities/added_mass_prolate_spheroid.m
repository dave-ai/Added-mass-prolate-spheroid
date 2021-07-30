%% Compute the added mass of a prolate spheroid
%
%
% author: Davide Grande
% date:   15/07/2021
%
% 
%
% Reference: 
% [1] This utility is part of the work published:
%   Open source simulation of underwater gliders, 
%   Davide Grande, L. Huang, C. A. Harris, P. Wu, G. Thomas1, E. Anderlini,
%   Oceans 2021, IEEE, San diego - Porto, September 2021.
%
%
% [2] THE COMPLETE EXPRESSIONS FOR "ADDED MASS" OF A RIGID BODY MOVING IN AN IDEAL FLUID
% Frederick H. Imlay, DAVID TAYLOR MODEL BASIN WASHINGTON DC
% Institution, Technical Report, July 1961.
% 
%

clear all 
close all 
clc


%% Parameters
b = 0.15; % minor semi-axis #1
c = 0.15; % minor semi-axis #2
a = 0.5;  % major semi-axis #3



mass_available = true; % Is the mass of the vehicle available? 
mom_of_inertia_available = true; % Is the moment of inertia of the vehicle available?
m = 56.34;   % [kg] -- mass
vol = 55.20; % [l]  -- volume 
I_yy = 2.62; % [kg m^2] -- moment of inertia around the y-axis 
I_zz = 2.62; % [kg m^2] -- moment of inertia around the z-axis 



%% Execution 
rho = m / vol*1000; % density of the vehicle 
e = sqrt(1-(b/a)^2); % eccentricity of the meridian elliptical section
e2 = e^2;

alpha_0 = 2*(1-e2)/e^3 * (1/2 * log((1+e)/(1-e)) - e);
beta_0 = 1 / e2 - (1-e2)/(2*e^3) * log((1+e)/(1-e));

%% Lam's factors
k1 = alpha_0 / (2-alpha_0);
k2 = beta_0 / (2-beta_0);
k_prime = e^4 * (beta_0 - alpha_0) / ( (2-e2) * (2*e2 - (2 - e^2)*(beta_0 - alpha_0) ) );


%% Computation of the added mass coefficients 
if (mass_available && mom_of_inertia_available)
    X_u_dot = -k1 * m;
    Y_v_dot = -k2 * m;
    Z_w_dot = -k2 * m;
    K_p_dot = 0;
    M_q_dot = -k_prime * I_yy;
    N_r_dot = -k_prime * I_zz;
    
else
    X_u_dot = -k1 * (4/3 * pi * rho * a * b^2);
    Y_v_dot = -k2 * (4/3 * pi * rho * a * b^2);
    Z_w_dot = -k2 * (4/3 * pi * rho * a * b^2);
    K_p_dot = 0;
    M_q_dot = -k_prime * 4 / 5 * pi * rho * a * b^2 * (a^2 + b^2);
    N_r_dot = -k_prime * 4 / 5 * pi * rho * a * b^2 * (a^2 + c^2);
end



%% Added mass matrix
M_A = diag([X_u_dot, Y_v_dot, Z_w_dot, K_p_dot, M_q_dot, N_r_dot])



