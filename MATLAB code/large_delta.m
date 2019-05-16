% Parameters
beta = 1;
delta = 1e4;

tmax = 100; % Maximum value of time
tstep = tmax / 100000; % Time steps
tvals = 0 : tstep : tmax; % Time domain

% Asymptotic expansions for large delta, beta << 1
s_inner = tvals.^2;
s_middle_altered = (2 / delta) * (1 - cos(sqrt(delta) * tvals));
% s_middle = (2 / delta) * (1 - cos(sqrt(delta) * tvals));


% Numerical calculation
[s, sdot, sddot] = numerical_solution(tvals, beta, delta);

B = - sqrt(delta) - atan((3 + beta) / (2 * sqrt(delta)));
A = -2 / cos(B + sqrt(delta));


sol = 2/delta + (A/delta) * cos(sqrt(delta) * sqrt(1 + 2 * tvals) + B) ...
    .* (1 + 2 * tvals).^(-(3 + beta) / 4);

B = -sqrt(delta);
A = -2;
simple_sol = 2/delta + (A/delta) * cos(sqrt(delta) * sqrt(1 + 2 * tvals) + B) ...
    .* (1 + 2 * tvals).^(-(3 + beta) / 4);
% Creating plots
close all

s_asy = (2/delta) * (1 - cos(sqrt(delta) * (sqrt(1 + 2 * tvals) - 1)) ...
.* (1 + 2 * tvals).^(-(3 + beta) / 4)); 

figure(1);
hold on
plot(tvals, s_inner);
plot(tvals, sol);
legend("Numerical", "Asymptotic");

figure(2);
plot(tvals, sdot);

