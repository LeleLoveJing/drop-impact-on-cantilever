%% 
% Comparison of asymptotic theory of an axisymmetric droplet impacting onto
% a cantilever with experiment.
%

epsilon = 0.1; % Error in non-dimensional value
rho = 789; % Density of fluid
M = 2.254e-5; % Mass of cantilever
U = 1.01; % Initial velocity of the droplet
R = 0.985e-3; % Initial radius of the droplet

% Reads the experimental data from the csv file
experiment_data = dlmread("RE_1328_displacement.csv");
experiment_t = experiment_data(:, 1); % Impact times
% Displacement expressed in terms of metres rather than millimetres
experiment_s = (experiment_data(:, 2) - experiment_data(1, 2)) * 10^-3;

% Theoretical predictions 
theory_t = 0:1e-5:max(experiment_t);
theory_s = (8 * sqrt(3) / 5) * (rho * U^(5/2) * R^(3/2) / M) * theory_t.^(5/2);

% Plotting
close all;
figure;
hold on;
scatter(experiment_t, experiment_s, 'filled');
plot(theory_t, theory_s);
xlabel("$t^*$/s", 'Interpreter', 'latex');
ylabel("$s^*(t^*)$/m", 'Interpreter', 'latex');
legend("Experimental data", "Asymptotic theory", 'Location','northwest');
title("Cantilever displacement");

figure;
hold on;
scatter(log(experiment_t), log(experiment_s), 'filled');
plot(log(theory_t), log(theory_s));
