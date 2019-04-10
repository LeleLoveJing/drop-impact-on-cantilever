%% Droplet impact on a cantilever: small_damping.m
% 
% This script deals with the distinguised limit where there is a small
% amount of damping, characterised by beta, which is larger than the square
% root of the mass ratio alpha but still much smaller than 1. It plots the
% first order correction to the cantilever displacement s_0(t) from the
% case where there is no damping.

%%
% Solution for cantilever displacement and its derivatives
beta = 1e-2; % Damping constant

s0 = @(t) 0.5 + t - 0.5 * sqrt(1 + 4 * t); % Leading order solution
s1 = @(t) ((1 + 4 * t) / 12 - ...
    (1 + 6 * t .* (1 + t)) ./ (12 * sqrt(1 + 4 * t))); % First order correction
s = @(t) s0(t) + beta * s1(t);

s = @(t) 0.5 + t - 0.5 * sqrt(1 + 4 * t) ...
    - (beta / 12) * ((1 + 6 * t + 6 * t.^2) ./ sqrt(1 + 4 * t) - 1 - 4 * t);
sdot = @(t) 1 - 1 ./ sqrt(1 + 4 * t) ...
    - (beta / 3) * ((1 + 6 * t + 9 * t.^2) / (1 + 4 * t).^(3/2) - 1);
sddot = @(t) 2 ./ (1 + 4 * t).^(3/2) ...
    - 2 * beta * (1 + 3 * t) ./ (1 + 4 * t).^(5/2);

%%
% Creating plots for a specific time array using the function
% plotting_cantilever
tvals = 0.01:0.01:10;
plotting_cantilever(tvals, s(tvals), sdot(tvals), sddot(tvals));
