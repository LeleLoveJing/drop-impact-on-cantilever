%% Droplet impact on a cantilever: small_damping_and_stiffness.m
% 
% This script deals with the distinguised limit where there is a small
% amount of damping, characterised by beta and a small amount of stiffness,
% characterised by delta, which is larger than the square root of the mass 
% ratio alpha but still much smaller than 1. It plots the first order 
% correction to the cantilever displacement s_0(t) to the case where there 
% is no damping or stiffness. 

%%
% Solution for cantilever displacement and its derivatives
beta = 1e-2;
delta = 1e-3;

% Homogeneous solution with no damping or stiffness
s0 = @(t) 0.5 + t - 0.5 * sqrt(1 + 4 * t); 
sdot0 = @(t) 1 - 1 ./ sqrt(1 + 4 * t);
sddot0 = @(t) 2 * (1 + 4 * t).^(-3/2);

s = @(t) s0(t) ...
    - (beta / 12) * ((1 + 6 * t + 6 * t.^2) ./ sqrt(1 + 4 * t) - 1 - 4 * t) ...
    - (delta / 120) ...
        * ((1 + 10 * t + 30 * t.^2 + 20 * t.^3) ./ sqrt(1 + 4 * t) ...
            - (1 + 4 * t).^2);
sdot = @(t) sdot0(t) ...
    - (beta / 3) * ((1 + 6 * t + 9 * t.^2) / (1 + 4 * t).^(3/2) - 1) ...
    - (delta / 15) ...
        * ((1 + 10 * t + 30 * t.^2 + 25 * t.^3) ./ (1 + 4 * t).^(3/2) ...
            - (1 + 4 * t));
sddot = @(t) sddot0(t) ...
    - 2 * beta * (1 + 3 * t) ./ (1 + 4 * t).^(5/2) ...
    - (delta / 15) ...
        * ((4 + 40 * t + 135 * t.^2 + 150 * t.^3) ./ (1 + 4 * t).^(5/2) - 4);
 
%%
% Later time solution for t > (12 / delta)^(2/3)
slate = @(t) t - (12/delta)^(2/3) * sqrt(5 * (delta / 12)^(2/3) * t - 4);

scomp = @(t) s0(t) + slate(t) - t;
%%
% Messing
tdelta = (12 / delta)^(2/3)


%%
% Creating plots for a specific time array using the function
% plotting_cantilever
tvals1 = 0.01:0.01:1.1 * tdelta;
tvals2 = 0.9 * tdelta:0.01:3 * tdelta;
tvalstot = 0.01:0.01:3 * tdelta;

hold on
plot(tvalstot, s(tvalstot));
plot(tvalstot, slate(tvalstot));
plot(tvalstot, scomp(tvalstot));
hold off
legend("s1", "s2", "scomp");

% plotting_cantilever(tvals, s(tvals), sdot(tvals), sddot(tvals));