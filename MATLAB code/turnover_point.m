function [d, ddot, dddot] = turnover_point(tvals, s, sdot, sddot)
% Calculates the turnover point given the solution for s

    d = 2 * sqrt(tvals - s);
    ddot = (1 - sdot) ./ sqrt(tvals - s);
    dddot = - 0.5 * ((1 - sdot).^2 +  2 * (tvals - s) .* sddot) ...
        ./ (2 * (tvals - s).^1.5);
end