function [s, sdot, sddot] = numerical_solution(tvals, beta, delta)
%NUMERICAL_SOLUTION Gives the numerical solution to the cantilever
%displacement for any order 1 damping beta and stiffness delta.
%   This function solves the governing ODE for the cantilever displacement
%   for any general beta and delta, given that they are order 1. It does
%   this by expressing the ODE as a system of first order ODEs and solving
%   using ode45. The input is tvals, the array for time which the ODE
%   should be solved over, and the physical parameters beta and delta. It 
%   outputs an array for the cantilever displacement s and its first and
%   second derivatives sdot and sddot.


% Function for the second derivative in terms of s and its derivative sdot
s_2nd_deriv = @(t, s, sdot) ...
    (2 * (1 - sdot).^2 - beta * sdot - delta * s) ...
        ./ (1 + 2 * (t - s));
 
% Function for solving the second order ODE as a system of first order
% ODEs, where s1 = s(t), s2 = s'(t), s_arr = [s1, s2] and ode_fun(t, s_arr)
% = [s1'(t), s2'(t)].
ode_fun = @(t, s_arr) ...
    [s_arr(2); s_2nd_deriv(t, s_arr(1), s_arr(2))];

% Solves the ODE over the time domain with initial conditions s(0) = s'(0)
% = 0.
[t, s_arr] = ode45(ode_fun, tvals, [0, 0]);

% Defines arrays for s and its derivatives
s = s_arr(:, 1);
sdot = s_arr(:, 2);
sddot = s_2nd_deriv(t, s, sdot);
end
