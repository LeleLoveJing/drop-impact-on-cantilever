%% Droplet impact on a cantilever: mass_only.m
% 
% This script deals with the distinguished limit where the mass term of the
% ODE for the cantilever dominates over the damping and spring term. In
% this case we have an analytic solution for s_0(t) for all time. We drop
% the 0 subscripts.

close all

%%
% Solution for s_0(t) and other quantities

% Cantilever displacement and its time derivatives
s = @(t) 0.5 + t - 0.5 * sqrt(1 + 4 * t); 
sdot = @(t) 1 - 1 ./ sqrt(1 + 4 * t);
sddot = @(t) 2 * (1 + 4 * t).^(-3/2);

% Turnover point and its time derivative 
d = @(t) 2 * sqrt(t - s(t));
ddot = @(t) (1 - sdot(t)) ./ sqrt(t - s(t));

% Jet thickness
hJ = @(t) (pi / 4) * (t - s(t)).^(3/2);


%%
% Specific parameters
tvals = 0.01:0.01:10;

%%
% Plotting cantilever displacement
figure(1);
plot(tvals, s(tvals));

figure(2);
plot(tvals, sdot(tvals));

figure(3);
plot(tvals, sddot(tvals));

%%
% Outer region solution

% Square root defined using the appropriate branch cut
branch_sqrt = @(zeta, t) abs(zeta - d(t)) .* abs(zeta + d(t)) ...
    .* exp(1i * (angle(zeta - d(t)) + angle(zeta + d(t))) / 2);

% Complex potential and its derivative
W = @(zeta, t) 1i * sdot(t) .* zeta + 1i * (1 - sdot(t)) .* branch_sqrt(zeta, t);
dW = @(zeta, t) 1i * sdot(t) + 1i * (1 - sdot(t)) .* zeta ./ branch_sqrt(zeta, t);

% Components of velocity in outer region
outer_u = @(X, Z, t) real(dW(X + 1i * Z, t));
outer_v = @(X, Z, t) -imag(dW(X + 1i * Z, t));

% Plots the velocity field
% x_axis = -6 : 0.1 : 6;
% z_axis = 0.01 : 0.01 : 1;
% [Xs, Zs] = meshgrid(x_axis, z_axis);
% figure(4);
% xlim([-6, 6]);
% ylim([0, 1]);
% for tval = tvals
%     us = outer_u(Xs, Zs, tval);
%     vs = outer_v(Xs, Zs, tval);
%    
%     figure(4);
%     scalefactor = 0.01;
%     quiver(Xs, Zs, scalefactor * us, scalefactor * vs, 'Autoscale', 'off');
% 
%     title(["t = ", num2str(tval), ", s(t) = ", num2str(s(tval))]);
%     pause(.1);
% end

% Free surface shape
H = @(X, t) 0.5 * abs(X) .* sqrt(X.^2 - 4 * (t - s(t))) - s(t);
tval = 0.1;
Xs = d(tval) * (1.0 : 1e-6 : 2);
figure(5);
pbaspect([1 1 1])
hold on
plot(Xs, H(Xs, tval));
plot(-Xs, H(-Xs, tval));
hold off
axis equal

% Pressure on the cantilever
figure(6);
P = @(X, t) - sddot(t) .* sqrt(d(t).^2 - X.^2) ... 
    + (1 - sdot(t)) .* d(t) .* ddot(t) ./ sqrt(d(t).^2 - X.^2);
tval = 0.1;
Xs = d(tval) * (-1 : 1e-3 : 1);
plot(Xs, P(Xs, tval));


%%
% Inner region solution

% Free surface shape 
surface_x = @(sigma, t) (hJ(t) / pi) * (sigma - log(sigma) -1);
inner_h = @(sigma, t) hJ(t) * (1 + 4 * sqrt(sigma) / pi);

sigmas = 0:1e-6:20;
tval = 0.1;
figure(7);
plot(surface_x(sigmas, tval), inner_h(sigmas, tval));



