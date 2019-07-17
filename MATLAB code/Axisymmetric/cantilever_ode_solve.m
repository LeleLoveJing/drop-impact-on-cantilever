
alpha = 30;
beta = 0.04;
gamma = 0.02;
epsilon = 1/10;

A = alpha / epsilon^3
B = beta / epsilon
C = epsilon * gamma

tvals = 0:1e-3:10;

s_2nd_deriv = @(t, s, sdot) ...
    (6 * sqrt(3) * sqrt(t - s) .* (1 - sdot.^2) - C * s ...
        - (B + 12 * sqrt(3) * sqrt(t - s)) .* sdot) ...
     ./ (A + 4 * sqrt(3) * (t - s).^(3/2));

ode_fun = @(t, s_arr) ...
    [s_arr(2); s_2nd_deriv(t, s_arr(1), s_arr(2))];

[t, s_arr] = ode45(ode_fun, tvals, [0, 0]);

% plot(tvals, s_arr(:, 1));

R = 0.985e-3;
U = 1.01;

plot(epsilon^2 * (R / U) * tvals, epsilon*2 * R * s_arr(:, 1) + 2.3e-3);