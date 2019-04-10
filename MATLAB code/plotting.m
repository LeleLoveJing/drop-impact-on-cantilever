function out = plotting(tvals, s, sdot, sddot)
%PLOTTING Creates figures of quantities of interest for droplet impact.
%   Given input of a time range tvals, along with arrays for the cantilever
%   displacement s and its first and second time derivative sdot and sddot,
%   this function creates plots of the different physical quantities that
%   we are interested in for the early time of droplet impact on a
%   cantilever.

close all % Closes all figures currently open

%%
% Important time dependent quantities

% Turnover point and its time derivatives
d = 2 * sqrt(tvals - s);
ddot = (1 - sdot) ./ sqrt(tvals - s);
dddot = - sddot ./ sqrt(tvals - s) ...
    - (1 - sdot).^2 ./ (2 * (tvals - s).^(3/2));

% Turnover point in the fixed plate case
d_fixed = 2 * sqrt(tvals);
ddot_fixed = 1 ./ sqrt(tvals);

% Jet thickness
hJ = (pi / 4) * (tvals - s).^(3/2);


%%
%  



end

