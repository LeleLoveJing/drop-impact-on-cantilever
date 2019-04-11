%% Droplet impact on a cantilever: small_damping.m
% 
% This script deals with the distinguised limit where there is a small
% amount of damping, characterised by beta, which is larger than the square
% root of the mass ratio alpha but still much smaller than 1. It plots the
% first order correction to the cantilever displacement s_0(t) from the
% case where there is no damping.

%%
% Solution for cantilever displacement and its derivatives
beta = 1e-1; % Damping constant

% Homogeneous solution with no damping
s0 = @(t) 0.5 + t - 0.5 * sqrt(1 + 4 * t); 
sdot0 = @(t) 1 - 1 ./ sqrt(1 + 4 * t);
sddot0 = @(t) 2 * (1 + 4 * t).^(-3/2);

% Solution with first order correction due to damping
s = @(t) s0(t) ...
    - (beta / 12) * ((1 + 6 * t + 6 * t.^2) ./ sqrt(1 + 4 * t) - 1 - 4 * t);
sdot = @(t) sdot0(t) ...
    - (beta / 3) * ((1 + 6 * t + 9 * t.^2) / (1 + 4 * t).^(3/2) - 1);
sddot = @(t) sddot0(t) ...
    - 2 * beta * (1 + 3 * t) ./ (1 + 4 * t).^(5/2);

%%
% Creating plots for a specific time array using the function
% plotting_cantilever
tvals = 0.01:1000:1e5;
close all
figure;
hold on
plot(tvals, s(tvals))
plot(tvals, s0(tvals));
plot(tvals, -beta * tvals.^1.5 / 4);
plot(tvals, tvals);
hold off
legend("Damped", "Mass only", "Limit", "Other limit");

% plotting_cantilever(tvals, s(tvals), sdot(tvals), sddot(tvals));
