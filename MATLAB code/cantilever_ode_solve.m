%% Droplet impact on a cantilever: cantilever_ode_solve.m
% 
% Main script for solving the ODE for the leading order displacement of the
% cantilever. Variables are expressed in dimensionless coordinates. We are
% solving the 2nd order non-linear ODE:
% (alpha / epsilson^2) s''(t) + beta s'(t) + epsilon^2 gamma s(t) = F(t)
% where F(t) is the hydrodynamic force, to leading order given by
% F_0(t) = 2 (-s_0''(t) (t - s_0(t)) + 2(1 - s_0(t))^2).
% 
% In this script, we solve the ODE using a specific epsilon, alpha, beta
% and gamma to leading order for s_0(t). From here on, we will use s(t) to
% mean s_0(t). 

%%
% Defining parameters
epsilon = 10^-4; % Small parameter for time
alpha = epsilon^2 ; % Ratio of mass of cantilever with mass of droplet
beta = 1; % Dimensionless damping constant
gamma = 1 / epsilon; % Dimensionless spring constant

total_time = 100; % Total dimensionless time to solve ODE for 

%%
% Solving the ODE: As this is a 2nd order ODE, we define the problem as a
% system of ODEs, with s1 = s(t) and s2 = s'(t). 
ode_fun = @(t, s_arr) ...
    [s_arr(2); ...
    (2 * (1 - s_arr(2))^2 - beta * s_arr(2) - epsilon^2 * gamma * s_arr(1)) ...
        / (alpha / epsilon^2 + 2 * (t - s_arr(1)))];
[t, s_arr] = ode45(ode_fun, [0, total_time], [0, 0]);

s = s_arr(:, 1); % Solution for s
sdot = s_arr(:, 2); % Time derivative of s

%%
% Now the ODE for s(t) is solved, can calculate other quantities
d = 2 * sqrt(t - s); % Turnover point d(t)

m = 2 * (1 - sdot) .* (t - s); % Equal to (d/dt)((t - s)^2)

% Second derivative of s
sddot = (2 * pi * (1 - sdot).^2 - beta * sdot - epsilon^2 * gamma * s) ...
    ./ (alpha^2 / epsilon^2 + 2 * pi * (t - s));

% Dimensionless force on the plate
F = 2 * pi * (- sddot .* (t - s) + (1 - sdot).^2);

% Energy in outer-outer region (dimensionless)
E_o_o = (1 - epsilon^2 * m);

% Energy in the outer region
E_o = epsilon^2 * sdot .* m;

% Energy in the jet region
E_J = zeros(length(t), 1);
for i = 1:length(t)
    % Numerically integrate with respect to time
    E_J(i) = trapz(t(1:i), (1 - sdot(1:i)).^3, 1);
end
E_J = 2 * epsilon^2 * E_J;


% Total work done
work_done = zeros(length(t), 1);
for i = 1:length(t)
    % Numerically integrate with respect to time
    work_done(i) = trapz(t(1:i), sdot(1:i) .* F(1:i), 1);
end

work_done = 2 * epsilon^2 * work_done;

% Total energy dissipated
energy_diss = zeros(length(t), 1);
for i = 1:length(t)
    % Numerically integrate with respect to time
    energy_diss(i) = trapz(t(1:i), sdot(1:i).^2, 1);
end

energy_diss = 2 * beta * epsilon^2 * energy_diss;

% Kinetic energy of the cantilever
E_K_cant = alpha * sdot.^2;

% Potential energy of the cantilever
E_P_cant = epsilon^4 * gamma * s.^2;

%%
% Plotting desired funtions
close all

% Cantilever displacement s(t)
figure;
plot(t, s)
title("s(t)");
xlabel("t")
ylabel("s(t)");
hold on
% plot(t, sqrt(epsilon^2 * gamma / 4) * t.^2);
hold off

% Derivative of s(t)
figure;
plot(t, sdot);
title("s'(t)");
xlabel("t");
ylabel("s'(t)");

% Second derivative of s(t)
figure;
plot(t, sddot);
title("s''(t)");
xlabel("t");
ylabel("s''(t)");

% t - s(t)
figure;
plot(t, t - s);
title("q(t)");
xlabel("t");
ylabel("t - s(t)");

% 1 - sdot 
figure;
plot(t, 1 - sdot);
title("q'(t)");
xlabel("t");
ylabel("1 - s'(t)");

% Turnover point d(t)
figure;
hold on
plot(t, d);
plot(t, 2 * sqrt(t));
hold off
title("d(t)");
xlabel("t");
ylabel("d(t)");
legend("Cantilever", "Fixed plate");

% m
figure;
plot(t, m);
title("(d / dt) ((t - s)^2)")
xlabel("t");
ylabel("m(t)");

% Energies
figure;
title("Energies");
hold on
plot(t, (E_o_o - 1) / epsilon^2);
plot(t, E_o / epsilon^2);
plot(t, E_J / epsilon^2);
plot(t, work_done / epsilon^2);
plot(t, - energy_diss / epsilon^2);
xlabel("t");
ylabel("Dimensionless energy / epsilon^2");
hold off
legend("E_{o-o} - 1", "E_o", "2 E_J", "Work done", "Energy dissipated");

% Jet energy
figure;
title("Kinetic energy of jet");
hold on
plot(t, E_J / epsilon^2);
plot(t, t);
hold off
xlabel("t");
ylabel("Dimensionless energy / epsilon^2");
legend("Cantilever", "Fixed plate");


% Energy balance
LHS = E_J - E_o - E_o_o + E_K_cant + E_P_cant + energy_diss;
RHS = 2 * epsilon^2 * m - 1;
figure;
title("Balance");
hold on
plot(t, LHS);
plot(t, RHS);
hold off
xlabel("t");
ylabel("Energy");
legend("LHS", "RHS");
