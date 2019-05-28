
close all


% Physical parameters
epsilon = 1e-3;
beta = 1;
delta = 1;

tvals = 0.01:0.01:10; % Time values

[s, sdot, sddot] = numerical_solution(tvals, beta, delta);
tvals = tvals';
[d, ddot, dddot] = turnover_point(tvals, s, sdot, sddot);

% Outer region pressure on cantilever
figure(1);
Xs = 1.1 * d(end) * (-1 : 1e-3 : 1);
Ps = - sddot(end) * sqrt(d(end)^2 - Xs.^2) ...
    + (1 - sdot(end)) * d(end) * ddot(end) ./ sqrt(d(end)^2 - Xs.^2);
plot(Xs, Ps,'LineWidth',2);
grid on;
xlabel('$X$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$P_0(X, 0, t)$', 'Interpreter', 'latex', 'FontSize', 16);
titletext = sprintf('Pressure on cantilever in outer region for $d_0(t)$ = %g', d(end));
title(titletext, 'Interpreter', 'latex', 'FontSize', 14);
print('figures/outer_pressure.png', '-dpng');

% Small beta and delta solution
[s_small, sdot_small, sddot_small] = asymptotic_solution(tvals, 0, 0);
figure(2);
plot(tvals, s_small, 'LineWidth',2);
grid on;
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$s_0(t)$', 'Interpreter', 'latex', 'FontSize', 16);
titletext = 'Displacement of cantilever for $\beta, \delta \ll 1$';
title(titletext, 'Interpreter', 'latex', 'FontSize', 14);
print('figures/cantilever_small_beta_and_delta.png', '-dpng');

% Beta and delta = 1
beta = 1;
delta = 1;

tvals = 0.01:0.01:100; % Time values

[s, sdot, sddot] = numerical_solution(tvals, beta, delta);



figure(3);
plot(tvals, s, 'LineWidth', 2);
grid on;
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$s_0(t)$', 'Interpreter', 'latex', 'FontSize', 16);
titletext = sprintf('Displacement of cantilever for $\\beta$ = %g and $\\delta$ = %g', beta, delta);
title(titletext, 'Interpreter', 'latex', 'FontSize', 14);
print('figures/cantilever_beta_delta_1.png', '-dpng');