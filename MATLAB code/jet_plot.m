% Parameters
beta = 1;
delta = 1e4;

epsilon = 1e-3;

tmax = 20; % Maximum value of time
tstep = tmax / 100000; % Time steps
tvals = 0 : tstep : tmax; % Time domain


% Numerical calculation
[s, sdot, sddot] = numerical_solution(tvals, beta, delta);
tvals = tvals';
[d, ddot, dddot] = turnover_point(tvals, s, sdot, sddot);

jet_xs = zeros(length(tvals), 1);
jet_hs = zeros(length(tvals), 1);



for idx =  1:length(tvals)
    tau = tvals(idx);
    jet_xs(idx) = 2 * ddot(idx) * (tmax - tau) + d(idx);
    jet_hs(idx) = (pi/4) * (ddot(idx) * (tau - s(idx))^1.5)...
        ./ (ddot(idx) - 2 * dddot(idx) * (tmax - tau));

end

figure;
plot(tvals, s);
xlabel("t");
ylabel("s(t)");
title("Cantilever displacement for delta >> 1");

figure;
plot(scaled_jet_xs, scaled_jet_hs);
xlabel("x");
ylabel("h");
title("Jet for delta >> 1 case");



d(end)