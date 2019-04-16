%% Droplet impact on a cantilever: general_damping_and_stiffness.m
% 
% This script deals with general values of damping and stiffness,
% characterised by beta and delta. There is no analytical solution in this
% case, so it uses ode45 to solve the differential equation numerically.
% This solution is only valid if beta and delta are O(1). 

%%
% Solution for cantilever displacement and its derivatives
beta = 1e-1; % Damping constant
delta = 1e-1; % Stiffness constant


