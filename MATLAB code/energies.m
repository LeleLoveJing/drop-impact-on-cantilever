function [outer_outer_energy, outer_energy, jet_energy, work_done, ...
    energy_diss] = energies(tvals, s, sdot, sddot, epsilon, beta)
%ENERGIES Calculates the different energies given a specific solution
%   This function takes input of a specific solution to the cantilever
%   problem, along with the governing physical parameters. It then
%   calculates the energies of interest for the problem. Specifically, it
%   calculates the kinetic energy of the droplet in the outer-outer region,
%   outer region and jet region, as well as the work done on the cantilever
%   and the energy dissipated. 

% Force on the plate, equal to (d^2/dt^2)[(t - s(t))^2]
force = 2 * ((1 - sdot).^2 - sddot .* (tvals - s));

% Energy in the outer-outer region
outer_outer_energy = 1 - 2 * epsilon^2 * (1 - sdot) .* (tvals - s);

% Energy in the outer region
outer_energy = - 2 * epsilon^2 * sdot .* (1 - sdot) .* (tvals - s);

% Energy in the jet region
jet_energy = zeros(length(tvals), 1);
for i = 1 : length(tvals)
    jet_energy(i) = 2 * epsilon^2 * trapz(tvals(1:i), (1 - sdot(1:i)).^3, 1);
end

% Work done on the cantilever
work_done = zeros(length(tvals), 1);
for i = 1 : length(tvals)
    work_done(i) = 2 * epsilon^2 * trapz(tvals(1:i), sdot(1:i) .* force(1:i), 1);
end

% Energy dissipated due to the damping
energy_diss = zeros(length(tvals), 1);
for i = 1 : length(tvals)
   energy_diss(i) = 2 * epsilon^2 * beta * trapz(tvals(1:i), sdot(1:i).^2, 1); 
end
