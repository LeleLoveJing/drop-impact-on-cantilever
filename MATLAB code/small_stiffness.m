%% Droplet impact on a cantilever: small_stiff.m
% 
% This script deals with the distinguised limit where there is a small
% amount of stiffness, characterised by delta, which is larger than the
% square root of the mass ratio alpha but still much smaller than 1. It
% plots the first order correction to the cantilever displacement s_0(t) to
% the case where there is no stiffness. 

%%
% Solution for cantilever displacement and its derivatives
delta = 1e-3; % Stiffness constant

% Homogeneous solution with no damping
s0 = @(t) 0.5 + t - 0.5 * sqrt(1 + 4 * t); 
sdot0 = @(t) 1 - 1 ./ sqrt(1 + 4 * t);
sddot0 = @(t) 2 * (1 + 4 * t).^(-3/2);


s = @(t) s0(t) - (delta / 120) ...
    * ((1 + 10 * t + 30 * t.^2 + 20 * t.^3) ./ sqrt(1 + 4 * t) ...
        - (1 + 4 * t).^2);
sdot = @(t) sdot0(t) - (delta / 15) ...
    * ((1 + 10 * t + 30 * t.^2 + 25 * t.^3) ./ (1 + 4 * t).^(3/2) ...
        - (1 + 4 * t));
sddot = @(t) sddot0(t) - (delta / 15) ...
    * ((4 + 40 * t + 135 * t.^2 + 150 * t.^3) ./ (1 + 4 * t).^(5/2) - 4);

%%
% Creating plots for a specific time array using the function
% plotting_cantilever
tvals = 0.01:0.01:100;

close all
figure;
hold on
plot(tvals, s(tvals))
plot(tvals, s0(tvals));
hold off
legend("Stiffness", "Mass only");
% plotting_cantilever(tvals, s(tvals), sdot(tvals), sddot(tvals));

