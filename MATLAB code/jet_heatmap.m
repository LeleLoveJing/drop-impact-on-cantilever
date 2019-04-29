%% Drop impact on a cantilever: jet_heatmap.m
% This script calculates the energy of the jet for a large number of
% cantilever parameters and then creates a heatmap of the results to
% display the energies. It will assess the energy at some large value of
% time where the validity of the solution starts to break down, i.e. as
% we're leaving the Wagner region. 

tmax = 30; % Maximum value of time

tvals = 0 : 1e-3 : tmax; % Time domain
tvals_trans = tvals'; % Transposes as used sometimes

% no_grid_points = 10;
% beta_max = 1;
% delta_max = 10;
% betas = 0 : beta_max / no_grid_points : beta_max; % Values of beta
% deltas = 0 : delta_max / no_grid_points : delta_max; % Values of delta 

betas = 0 : 0.1 : 1;
deltas = 0 : 0.1 : 10;

epsilon = 1; % Needs fixing, just sets so we get the answer without 

jet_energies = zeros(length(betas), length(deltas)); % Array for energies

% Function for the final jet energy given a solution for sdot
final_jet_energy = @(sdot) 2 * trapz(tvals, (1 - sdot).^3, 1);

% Iterates over all the given values of beta and delta
for beta_idx = 1 : length(betas)
    for delta_idx = 1 : length(deltas)
        % Solve the cantilever ODE for the given beta and delta
        disp([beta_idx, delta_idx]);
        [s, sdot, sddot] = ...
            numerical_solution(tvals, betas(beta_idx), deltas(delta_idx));
        disp("Solved ODE");
        
        % Calculates and saves the jet energy
        jet_energies(beta_idx, delta_idx) ...
            = final_jet_energy(sdot);
        disp("Solved energy");
    end
end

% Plots the colour map
figure;
imagesc(betas, deltas, jet_energies);
hold on;
% contour(betas, deltas, jet_energies, 'LineColor', 'k');
set(gca,'YDir','normal')
colorbar
xlabel("Damping term, beta")
ylabel("Spring constant, delta");