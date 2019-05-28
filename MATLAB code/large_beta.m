
% Parameters
beta = 1e3;
delta = 1e2;

tmax = 100; % Maximum value of time
tstep = tmax / 1000; % Time steps
tvals = 0 : tstep : tmax; % Time domain

% Asymptotic expansions for large beta
s_inner = 2 * (tvals / beta + exp(-beta * tvals) / beta^2 - 1 / beta^2);
s_outer =  (2 / delta)  * (1 - exp(-delta * tvals / beta));
s_composite = (2 / delta) * (1 - exp(-delta * tvals / beta)) ...
    -(2 / beta^2) * (1 - exp(-beta * tvals));
% jet_energy_asymptotic = (1 - 6 / beta) * tvals;
jet_energy_asymptotic = tvals - 6 / delta;

% Numerical calculation
[s, sdot, sddot] = numerical_solution(tvals, beta, delta);
jet_energy_numerical = zeros(length(tvals), 1);
for i = 2 : length(tvals)
    jet_energy_numerical(i) = jet_energy_numerical(i-1) ...
        + trapz(tvals(i-1:i), (1 - sdot(i-1:i)).^3, 1);
end

% Creating plots
close all

figure(1);
plot(tvals, s_composite, 'LineWidth',2);
grid on;
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$s_0(t)$', 'Interpreter', 'latex', 'FontSize', 16);
titletext = sprintf('Displacement of cantilever for $\\beta$ = %g and $\\delta$ = %g', beta, delta);
title(titletext, 'Interpreter', 'latex', 'FontSize', 14);
print('figures/large_beta.png', '-dpng');

% figure(2);
% plot(tvals, (s' - s_composite) );
% xlabel("t");
% ylabel("Error");
% title("Error for s(t)");
% 
% figure(3);
% hold on
% plot(tvals, jet_energy_numerical);
% plot(tvals, jet_energy_asymptotic);
% xlabel("t");
% ylabel("Jet energy");
% legend("Numerical", "Asymptotic");
% 
% figure(4);
% plot(tvals, jet_energy_numerical' - jet_energy_asymptotic);
% xlabel("t");
% ylabel("Error");
% title("Error for jet energy");