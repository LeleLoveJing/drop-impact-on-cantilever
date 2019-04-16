%% Drop impact on a cantilever: main.m
% This file is a general script for plotting and computing desired
% quantities on the fly. There's no real structure but results from here
% should translate themselves into separate scripts at some point

%% Validating asymptotics
% This section is just to validate that the asymptotic solutions are
% accurate for small delta and beta, by plotting the results against those
% found numerically. 

beta = 1e-3; % Damping coefficient
delta = 1e-3; % Stiffness coefficient

tvals = 0.01:0.001:100; % Time range

% Asymptotic solution
[s_asy, sdot_asy, sddot_asy] = asymptotic_solution(tvals, beta, delta);

% Numerical solution
[s_num, sdot_num, sddot_num] = numerical_solution(tvals, beta, delta);

% Homogeneous solution (no damping or stiffness)
[s_homo, sdot_homo, sddot_homo] = asymptotic_solution(tvals, 0, 0);

% Plotting s
figure(1);
hold on;
plot(tvals, s_asy);
plot(tvals, s_num);
plot(tvals, s_homo);
hold off;
legend("Asymptotic", "Numerical", "Homogeneous");