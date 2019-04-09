%% Droplet impact on a cantilever: mass_only.m
% 
% This script deals with the distinguished limit where the mass term of the
% ODE for the cantilever dominates over the damping and spring term. In
% this case we have an analytic solution for s_0(t) for all time. We drop
% the 0 subscripts.

close all % Closes all figures currently open

%%
% Solution for s_0(t) and other quantities

% Cantilever displacement and its time derivatives
s = @(t) 0.5 + t - 0.5 * sqrt(1 + 4 * t); 
sdot = @(t) 1 - 1 ./ sqrt(1 + 4 * t);
sddot = @(t) 2 * (1 + 4 * t).^(-3/2);

% Turnover point and its time derivative 
d = @(t) 2 * sqrt(t - s(t));
ddot = @(t) (1 - sdot(t)) ./ sqrt(t - s(t));
dddot = @(t) - sddot(t) ./ sqrt(t - s(t)) ...
    - (1 - sdot(t)).^2 ./ (2 * (t - s(t)).^(3/2));

% Turnover point in the case of a fixed plate
d_fixed = @(t) 2 * sqrt(t);
ddot_fixed = @(t) 1 ./ sqrt(t);

% Jet thickness
hJ = @(t) (pi / 4) * (t - s(t)).^(3/2);


%%
% Specific parameters
tvals = 0.01:0.01:10;

%%
% Plotting cantilever displacement
plotno = 1;
figure(plotno);
plot(tvals, s(tvals));
xlabel("t");
ylabel("s(t)");
title("Cantilever displacement, s(t)");

plotno = plotno + 1;
figure(plotno);
plot(tvals, sdot(tvals));
xlabel("t");
ylabel("s'(t)");
title("Time derivative of cantilever displacement, s'(t)");

plotno = plotno + 1;
figure(plotno);
plot(tvals, sddot(tvals));
xlabel("t");
ylabel("s''(t)");
title("Second time derivative of cantilever displacement, s''(t)");

%%
% Turnover point evolution
plotno = plotno + 1;
figure(plotno);
hold on
plot(tvals, d(tvals));
plot(tvals, 2 * sqrt(tvals));
xlabel("t");
ylabel("d(t)");
legend("Cantilever", "Fixed plate");
title("Location of turnover point, d(t)");

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
plotno = plotno + 1;
figure(plotno);
pbaspect([1 1 1])
hold on
plot(Xs, H(Xs, tval), 'b');
plot(-Xs, H(-Xs, tval), 'b');
hold off
xlim(1.1 * [-max(Xs), max(Xs)]);
ylim(1.1 * [0, max(H(Xs, tval))]);
daspect([1 1 1]);
xlabel("X");
ylabel("Z");
title(sprintf("Free surface H(X, t) in outer region at time t = %g", tval)); 

% Pressure on the cantilever
P = @(X, t) - sddot(t) .* sqrt(d(t).^2 - X.^2) ... 
    + (1 - sdot(t)) .* d(t) .* ddot(t) ./ sqrt(d(t).^2 - X.^2);
P_fixed = @(X, t) d_fixed(t) .* ddot_fixed(t) ./ sqrt(d_fixed(t).^2 - X.^2);

tval = 0.1;
Xs = d(tval) * (-1 : 1e-3 : 1);
Xs_fixed = d_fixed(tval) * (-1 : 1e-3 : 1);

plotno = plotno+1;
figure(plotno);
hold on
plot(Xs, P(Xs, tval));
plot(Xs, P_fixed(Xs_fixed, tval));
xlabel("X");
ylabel("P(X, 0, t)");
legend("Cantilever", "Fixed plate");
title(sprintf("Pressure on cantilever in outer region at time t = %g", tval));


%%
% Inner region solution

% Free surface shape 
surface_x = @(sigma, t) (hJ(t) / pi) * (sigma - log(sigma) -1);
inner_h = @(sigma, t) hJ(t) * (1 + 4 * sqrt(sigma) / pi);

sigmas = 0:1e-6:20;
tval = 0.1;
plotno = plotno + 1;
figure(plotno);
plot(surface_x(sigmas, tval), inner_h(sigmas, tval));
xlabel("$$\hat{x}$$", 'Interpreter','latex', 'FontSize', 18);
ylabel("$$\hat{z}$$", 'Interpreter','latex', 'FontSize', 18);
title(sprintf("Free surface shape in the inner region at time t = %g", tval));

% Pressure on the cantilever
cantilever_x = @(sigma, t) ...
    -(hJ(t) / pi) * (sigma + 4 * sqrt(sigma) + log(sigma) + 1);
inner_p = @(sigma, t) 2 * ddot(t)^2 * sqrt(sigma) ./ (1 + sqrt(sigma)).^2;

plotno = plotno + 1;
figure(plotno);
sigmas_1 = 0:1e-8:0.5;
hold on
plot(cantilever_x(sigmas_1, tval), inner_p(sigmas_1, tval), 'b');
sigmas_2 = 0.5:1e-4:20;
plot(cantilever_x(sigmas_2, tval), inner_p(sigmas_2, tval), 'b');
hold off
xlabel("$$\hat{x}$$", 'Interpreter','latex', 'FontSize', 18);
ylabel("$$\hat{p}_0(\hat{x}, 0, t)$$", 'Interpreter','latex', 'FontSize', 18);
title(sprintf("Pressure on the cantilever in the inner region at time t = %g", tval));


%%
% Jet region solution
jet_x = @(tau, t) 2 * ddot(tau) .* (t - tau) .* d(tau);
jet_u = @(tau) 2 * ddot(tau);
jet_h = @(tau, t) (pi / 4) * ddot(tau) .* (tau - s(tau)).^(3/2) ...
    ./ (ddot(tau) - 2 * dddot(tau) .* (t - tau));

% Free surface plot
tval = 0.1;
taus = 0:1e-6:tval;
plotno = plotno + 1;
figure(plotno);
plot(jet_x(taus, tval), jet_h(taus, tval));
xlim(1.1 * [0, max(jet_x(taus, tval))]);
ylim(2 * [0, max(jet_h(taus, tval))]);
daspect([1 1 1]);
xlabel("$$\bar{x}$$", 'Interpreter', 'latex', 'FontSize', 18);
ylabel("$$\bar{z}$$", 'Interpreter', 'latex', 'FontSize', 18);



