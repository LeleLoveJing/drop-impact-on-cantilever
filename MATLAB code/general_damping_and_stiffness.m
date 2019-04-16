%% Droplet impact on a cantilever: general_damping_and_stiffness.m
% 
% This script deals with general values of damping and stiffness,
% characterised by beta and delta. There is no analytical solution in this
% case, so it uses ode45 to solve the differential equation numerically.
% This solution is only valid if beta and delta are O(1). Recall the ODE
% that is being solved is
% s''(t) + beta * s'(t) + delta * s(t) = (d^2/dt^2)[(t - s(t))^2],
% which can be rearranged to give a formula for s''(t)
% s''(t) = [2(1 - s'(t))^2 - beta * s'(t) - delta * s(t)] ...
%              / (1 + 2(t - s(t))),
% The initial conditions are s(0) = s'(0) = 0.


%%
% Solution for cantilever displacement and its derivatives
beta = 1e-3; % Damping constant
delta = 1e-3; % Stiffness constant
tvals = 0.01:0.01:100; % Time domain

% Function for the second derivative in terms of s and its derivative sdot
s_2nd_deriv = @(t, s, sdot) ...
    (2 * (1 - sdot).^2 - beta * sdot - delta * s) ...
        ./ (1 + 2 * (t - s));
 
% Function for solving the second order ODE as a system of first order
% ODEs, where s1 = s(t), s2 = s'(t), s_arr = [s1, s2] and ode_fun(t, s_arr)
% = [s1'(t), s2'(t)].
ode_fun = @(t, s_arr) ...
    [s_arr(2); s_2nd_deriv(t, s_arr(1), s_arr(2))];

% Solves the ODE over the time domain
[t, s_arr] = ode45(ode_fun, tvals, [0, 0]);

% Defines arrays for s and its derivatives
s = s_arr(:, 1);
sdot = s_arr(:, 2);
sddot = s_2nd_deriv(t, s, sdot);


%%
% Plotting
plot(tvals, s);
