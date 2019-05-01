tmax = 2e2; % Maximum value of time
tvals = 0 : 1e-3 : tmax; % Time domain

beta = 1e2;
delta = 1e-3;

s_inner = 2 * (tvals / beta + exp(-beta * tvals) / beta^2 - 1 / beta^2);
% s_middle = (2 / beta) * (tvals - 1/beta);
s_outer =  (2 / delta)  * (1 - exp(-delta * tvals / beta));
s_composite = (2 / beta^2) * (exp(-beta * tvals) - 1) ...
    + (2 / delta) * (1 - exp(-delta * tvals / beta));

close all
figure;
hold on
[s, sdot, sddot] = numerical_solution(tvals, beta, delta);
plot(tvals, s);
plot(tvals, s_inner);
% plot(tvals, s_middle);
% plot(tvals, s_outer);
% plot(tvals, s_composite);
legend("Numerical", "Inner");


figure;
plot(tvals, s' - s_inner);

