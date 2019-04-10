function out = plotting_cantilever(tvals, s, sdot, sddot)
%PLOTTING_CANTILEVER Creates figures of quantities of interest for droplet impact.
%   Given input of a time range tvals, along with arrays for the cantilever
%   displacement s and its first and second time derivative sdot and sddot,
%   this function creates plots of the different physical quantities that
%   we are interested in for the early time of droplet impact on a
%   cantilever.

close all % Closes all figures currently open
plotno = 0; % Tracks how many plots there are


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
% Values of the time dependent functions at a specific time of our choosing
% in order to plot a slice of the physical quantities at those times

tstep = 12; % The timestep which we're plotting at
tval = tvals(tstep); % The value of time at the timestep we're plotting at
sval = s(tstep); % Value of cantilever displacement at tval
sdotval = sdot(tstep);
sddotval = sddot(tstep);
dval = d(tstep);
ddotval = ddot(tstep);
dddotval = dddot(tstep);
hJval = hJ(tstep);


%%
% Plotting cantilever displacement
plotno = plotno + 1;
figure(plotno);
plot(tvals, s);
xlabel("t");
ylabel("s(t)");
title("Cantilever displacement, s(t)");

plotno = plotno + 1;
figure(plotno);
plot(tvals, sdot);
xlabel("t");
ylabel("s'(t)");
title("Time derivative of cantilever displacement, s'(t)");

plotno = plotno + 1;
figure(plotno);
plot(tvals, sddot);
xlabel("t");
ylabel("s''(t)");
title("Second time derivative of cantilever displacement, s''(t)");

%%
% Turnover point evolution
plotno = plotno + 1;
figure(plotno);
hold on
plot(tvals, d);
plot(tvals, d_fixed);
xlabel("t");
ylabel("d(t)");
legend("Cantilever", "Fixed plate");
title("Location of turnover point, d(t)");

%%
% Outer region solution

% Velocity plotting: RETURN LATER

% % Square root defined using the appropriate branch cut
% branch_sqrt = @(zeta, t) abs(zeta - d(t)) .* abs(zeta + d(t)) ...
%     .* exp(1i * (angle(zeta - d(t)) + angle(zeta + d(t))) / 2);
% 
% % Complex potential and its derivative
% W = @(zeta, t) 1i * sdot(t) .* zeta + 1i * (1 - sdot(t)) .* branch_sqrt(zeta, t);
% dW = @(zeta, t) 1i * sdot(t) + 1i * (1 - sdot(t)) .* zeta ./ branch_sqrt(zeta, t);
% 
% % Components of velocity in outer region
% outer_u = @(X, Z, t) real(dW(X + 1i * Z, t));
% outer_v = @(X, Z, t) -imag(dW(X + 1i * Z, t));
%
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
H = @(X, sval, tval) 0.5 * abs(X) .* sqrt(X.^2 - 4 * (tval - sval)) - sval;

Xs = dval * (1.0 : 1e-6 : 2);
plotno = plotno + 1;
figure(plotno);
hold on
plot(Xs, H(Xs, sval, tval), 'b');
plot(-Xs, H(-Xs, sval, tval), 'b');
hold off
xlim(1.1 * [-max(Xs), max(Xs)]);
ylim(1.1 * [0, max(H(Xs, sval, tval))]);
daspect([1 1 1]);
xlabel("X");
ylabel("Z");
title(sprintf("Free surface H(X, t) in outer region at time t = %g", tval)); 

% Pressure on the cantilever
P = @(X, dval, ddotval, sdotval, sddotval, tval) - sddotval .* sqrt(dval.^2 - X.^2) ... 
    + (1 - sdotval) .* dval .* ddotval ./ sqrt(dval.^2 - X.^2);

P_fixed = @(X, t) d_fixed(tstep) .* ddot_fixed(tstep) ./ sqrt(d_fixed(tstep).^2 - X.^2);


Xs = dval * (-1 : 1e-3 : 1);
Xs_fixed = d_fixed(tstep) * (-1 : 1e-3 : 1);

plotno = plotno+1;
figure(plotno);
hold on
plot(Xs, P(Xs, dval, ddotval, sdotval, sddotval, tval));
plot(Xs, P_fixed(Xs_fixed, tval));
xlabel("X");
ylabel("P(X, 0, t)");
legend("Cantilever", "Fixed plate");
title(sprintf("Pressure on cantilever in outer region at time t = %g", tval));

%%
% Inner region solution

% Free surface shape 
surface_x = @(sigma, hJval) (hJval / pi) * (sigma - log(sigma) -1);
inner_h = @(sigma, hJval) hJval * (1 + 4 * sqrt(sigma) / pi);

sigmas = 0:1e-6:20;
plotno = plotno + 1;
figure(plotno);
plot(surface_x(sigmas, hJval), inner_h(sigmas, hJval));
xlabel("$$\hat{x}$$", 'Interpreter','latex', 'FontSize', 18);
ylabel("$$\hat{z}$$", 'Interpreter','latex', 'FontSize', 18);
title(sprintf("Free surface shape in the inner region at time t = %g", tval));

% Pressure on the cantilever
cantilever_x = @(sigma, hJval) ...
    -(hJval / pi) * (sigma + 4 * sqrt(sigma) + log(sigma) + 1);
inner_p = @(sigma, ddotval) 2 * ddotval^2 * sqrt(sigma) ./ (1 + sqrt(sigma)).^2;

plotno = plotno + 1;
figure(plotno);
sigmas_1 = 0:1e-7:1e-4;
hold on
plot(cantilever_x(sigmas_1, hJval), inner_p(sigmas_1, ddotval), 'b');
sigmas_2 = 1e-4:1e-4:20;
plot(cantilever_x(sigmas_2, hJval), inner_p(sigmas_2, ddotval), 'b');
hold off
xlabel("$$\hat{x}$$", 'Interpreter','latex', 'FontSize', 18);
ylabel("$$\hat{p}_0(\hat{x}, 0, t)$$", 'Interpreter','latex', 'FontSize', 18);
title(sprintf("Pressure on the cantilever in the inner region at time t = %g", tval));


%%
% Jet region solution
tau_range = tvals(1:tstep);
d_range = d(1:tstep);
ddot_range = ddot(1:tstep);
dddot_range = dddot(1:tstep);
s_range = s(1:tstep);

jet_x = 2 * ddot_range .* (tval - tau_range) .* d_range;
jet_u = 2 * ddot_range;
jet_h = (pi / 4) * ddot_range .* (tau_range - s_range).^(3/2) ...
    ./ (ddot_range - 2 * dddot_range .* (tval - tau_range));

plotno = plotno + 1;
figure(plotno);
plot(jet_x, jet_h);
xlim(1.1 * [0, max(jet_x)]);
ylim(2 * [0, max(jet_h)]);
daspect([1 1 1]);
xlabel("$$\bar{x}$$", 'Interpreter', 'latex', 'FontSize', 18);
ylabel("$$\bar{z}$$", 'Interpreter', 'latex', 'FontSize', 18);
title(sprintf("Free surface of jet at time t = %g", tval));


end

