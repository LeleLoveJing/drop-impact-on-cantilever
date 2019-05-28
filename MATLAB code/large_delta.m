% Parameters
beta = 0;
delta = 1e4;

tmax = 15; % Maximum value of time
tstep = tmax / 100000; % Time steps
tvals = 0 : tstep : tmax; % Time domain


% Numerical calculation
[s, sdot, sddot] = numerical_solution(tvals, beta, delta);

B = - sqrt(delta) - atan((3 + beta) / (2 * sqrt(delta)));
A = -2 / cos(B + sqrt(delta));


s_asy = 2/delta + (A/delta) * cos(sqrt(delta) * sqrt(1 + 2 * tvals) + B) ...
    .* (1 + 2 * tvals).^(-(3 + beta) / 4);

% Creating plots
close all

figure(1);
plot(tvals, s_asy, 'LineWidth',1);
grid on;
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$s_0(t)$', 'Interpreter', 'latex', 'FontSize', 16);
titletext = sprintf('Displacement of cantilever for $\\beta$ = %g and $\\delta$ = %g', beta, delta);
title(titletext, 'Interpreter', 'latex', 'FontSize', 11.5);
print('figures/large_delta.png', '-dpng');



