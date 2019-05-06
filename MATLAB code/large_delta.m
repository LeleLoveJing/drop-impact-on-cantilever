% Parameters
beta = 0;
delta = 1e4;

tmax = 1000 / sqrt(delta) ; % Maximum value of time
tstep = tmax / 1000; % Time steps
tvals = 0 : tstep : tmax; % Time domain

% Asymptotic expansions for large delta, beta << 1
s_inner = tvals.^2;
s_middle = (2 / delta) * (1 - cos(sqrt(delta) * tvals));

% Numerical calculation
[s, sdot, sddot] = numerical_solution(tvals, beta, delta);


% Creating plots
close all

figure(1);
hold on
plot(tvals, s);
% plot(tvals, s_middle);
hold off
xlabel("t");
ylabel("s(t)");
legend("Numerical", "Asymptotic");