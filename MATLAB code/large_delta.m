% Parameters
beta = 1;
delta = 1e4;

tmax = 50; % Maximum value of time
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
hold on
plot(tvals, s);
plot(tvals, s_asy);
legend("Numerical", "Asymptotic");
xlabel("t");
ylabel("s(t)");



