%% Droplet impact on a cantilever: mass_only.m
% 
% This script deals with the distinguished limit where the mass term of the
% ODE for the cantilever dominates over the damping and spring term. In
% this case we have an analytic solution for s_0(t) for all time. We drop
% the 0 subscripts.

close all % Closes all figures currently open

%%
% Solution for cantilever displacement and its derivatives
s = @(t) 0.5 + t - 0.5 * sqrt(1 + 4 * t); 
sdot = @(t) 1 - 1 ./ sqrt(1 + 4 * t);
sddot = @(t) 2 * (1 + 4 * t).^(-3/2);

%%
% Creating plots for a specific time array using the function "plotting"
tvals = 0.01:0.01:10;
plotting(tvals, s(tvals), sdot(tvals), sddot(tvals));

