%% Drop impact on a cantilever: main.m
% This file is a general script for plotting and computing desired
% quantities on the fly. There's no real structure but results from here
% should translate themselves into separate scripts at some point

%% Validating asymptotics
% This section is just to validate that the asymptotic solutions are
% accurate for small delta and beta, by plotting the results against those
% found numerically. 

% beta = 1e-3; % Damping coefficient
% delta = 1e-3; % Stiffness coefficient
% 
% tvals = 0.01:0.001:100; % Time range
% 
% % Asymptotic solution
% [s_asy, sdot_asy, sddot_asy] = asymptotic_solution(tvals, beta, delta);
% 
% % Numerical solution
% [s_num, sdot_num, sddot_num] = numerical_solution(tvals, beta, delta);
% 
% % Homogeneous solution (no damping or stiffness)
% [s_homo, sdot_homo, sddot_homo] = asymptotic_solution(tvals, 0, 0);
% 
% % Plotting s
% figure(1);
% hold on;
% plot(tvals, s_asy);
% plot(tvals, s_num);
% plot(tvals, s_homo);
% hold off;
% legend("Asymptotic", "Numerical", "Homogeneous");

%% Energy calculations
% This section is to calculate the energies of a specific problem, plot
% them, and see if they are as we expect. 

% Physical parameters
epsilon = 1e-4;
beta = 1;
delta = 1;

tvals = 0.01:0.01:4; % Time values

% Solves the problem numerically
[s, sdot, sddot] = numerical_solution(tvals, beta, delta);
figure(1)
plot(tvals, s)

tvals = tvals'; % Transposes for ease of addition

% Calculates the energies
[outer_outer_energy, outer_energy, jet_energy, work_done, energy_diss] ...
    = energies(tvals, s, sdot, sddot, epsilon, beta);

drop_energy = outer_outer_energy + outer_energy + jet_energy;
cantilever_energy = work_done - energy_diss;
total_energy = drop_energy + cantilever_energy;

figure(2);
hold on
plot(tvals, outer_outer_energy - 1);
plot(tvals, outer_energy);
plot(tvals, jet_energy);
plot(tvals, work_done);
plot(tvals, -energy_diss);
plot(tvals, total_energy  - 1);
hold off
legend("Outer-outer - 1", "Outer", "Jet", "Work done", "Energy dissipated", "Total - 1");
