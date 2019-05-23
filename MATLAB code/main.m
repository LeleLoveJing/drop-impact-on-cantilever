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
% 
close all

% Physical parameters
epsilon = 1e-3;
beta = 1;
delta = 10;

alpha = epsilon^2;
gamma = beta / epsilon^2;

tvals = 0.01:0.01:100; % Time values

% Solves the problem numerically
[s, sdot, sddot] = numerical_solution(tvals, beta, delta);
figure(1);
plot(tvals, s);
xlabel("t");
ylabel("s(t)");
title("Cantilever displacement s(t)");

tvals = tvals'; % Transposes for ease of addition

% Calculates the energies
[outer_outer_energy, outer_energy, jet_energy, cant_kinetic, ...
    cant_potential, work_done, energy_diss] ...
    = energies(tvals, s, sdot, sddot, epsilon, beta, delta);

total_drop_energy = outer_outer_energy + outer_energy + jet_energy;
total_cantilever_energy = cant_kinetic + cant_potential;
total_energy = total_drop_energy + total_cantilever_energy;

figure(2);
hold on;
plot(tvals, total_drop_energy - 1);
plot(tvals, work_done);
legend("Droplet energy loss", "Work done on cantilever");
xlabel("t");
ylabel("Energy");
title("Energy loss of droplet");


figure(3);
hold on;
plot(tvals, total_cantilever_energy);
plot(tvals, cant_kinetic);
plot(tvals, cant_potential);
plot(tvals, energy_diss);
legend("Total cantilever energy", "Kinetic energy", "Potential energy", "Dissipation");
xlabel("t");
ylabel("Energy");
title("Energy of cantilever");

figure(4);
hold on;
plot(tvals, total_cantilever_energy);
plot(tvals, energy_diss);
plot(tvals, total_drop_energy - 1);
xlabel("t");
ylabel("Energy");
legend("Total cantilever energy", "Energy dissipation", "Droplet energy loss");
title("Droplet-cantilever energy comparison");

figure(5);
hold on;
plot(tvals, jet_energy); 
plot(tvals, outer_outer_energy + outer_energy - 1);
xlabel("t");
ylabel("Energy");
legend("Jet energy");
legend("Jet energy", "Energy loss from bulk of droplet");
title("Energy of droplet");

%% Jet energy
% Plays around with the jet kinetic energy for different values of beta and
% delta 

% close all
% 
% % Physical parameters
% epsilon = 1e-3;
% 
% tvals = 0:0.01:30; % Time values
% 
% 
% figure(1);
% hold on
% title("s(t)");
% xlabel("t");
% ylabel("s(t)");
% legend show
% legend("Location", "best");
% 
% figure(2);
% hold on
% title("s'(t)");
% xlabel("t");
% ylabel("s'(t)");
% legend show
% legend("Location", "best");
% 
% figure(3);
% hold on
% title("Jet energy");
% xlabel("t");
% ylabel("E_J");
% legend show
% legend("Location", "best");
% 
% figure(4);
% hold on
% title("Work done");
% xlabel("t");
% ylabel("E_{work}");
% legend show
% legend("Location", "best");
% 
% for delta = 1
% for beta = 1
%     
%     plotname = sprintf("beta = %0.1e, delta = %0.1e", beta, delta);
% % Solves the problem numerically
%     [s, sdot, sddot] = numerical_solution(tvals, beta, delta);
%     
%     figure(1);
%     plot(tvals, s, 'DisplayName', plotname);
%     
%     figure(2);
%     plot(tvals, sdot, 'DisplayName', plotname);
%     
%     tvals = tvals'; % Transposes for ease of addition
% 
%     % Calculates the energies
%     [outer_outer_energy, outer_energy, jet_energy, work_done, energy_diss] ...
%         = energies(tvals, s, sdot, sddot, epsilon, beta);
%     
%     tvals = tvals';
%     
%     figure(3);
%     plot(tvals, (1/epsilon^2) * jet_energy, 'DisplayName', plotname);
%     
%     figure(4);
%     plot(tvals, (1/epsilon^2) * work_done, 'DisplayName', plotname);
% end
% end

