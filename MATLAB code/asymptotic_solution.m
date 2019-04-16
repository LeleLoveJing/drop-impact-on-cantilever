function [s, sdot, sddot] = asymptotic_solution(tvals, beta, delta)
%ASYMPTOTIC_SOLUTION Gives the asymptotic solution to the cantilever
%displacement for small damping beta and stiffness delta.
%   This function deals with the distinguised limit where there is a small
%   amount of damping, characterised by beta, and a small amount of
%   stiffness, characterised by delta, which are both larger than the
%   square root of the mass ratio, characterised by sqrt(alpha). The
%   solutions were derived using matched asymptotic expansions, so are in
%   the form of first order corrections to a homogeneous solution,
%   proportional to beta and delta. The input is tvals, an array for times
%   which want to be plotted, and then the physical parameters beta and
%   delta. It outputs an array for the cantilever displacement s and its
%   first and second derivatives sdot and sddot. 

% Functional homogeneous solution with no damping or stiffness
s0fun = @(t) 0.5 + t - 0.5 * sqrt(1 + 4 * t); 
sdot0fun = @(t) 1 - 1 ./ sqrt(1 + 4 * t);
sddot0fun = @(t) 2 * (1 + 4 * t).^(-3/2);

% Functional solution for solution with first order corrections for damping
% and stiffness
sfun = @(t) s0fun(t) ...
    - (beta / 12) * ((1 + 6 * t + 6 * t.^2) ./ sqrt(1 + 4 * t) - 1 - 4 * t) ...
    - (delta / 120) ...
        * ((1 + 10 * t + 30 * t.^2 + 20 * t.^3) ./ sqrt(1 + 4 * t) ...
            - (1 + 4 * t).^2);
sdotfun = @(t) sdot0fun(t) ...
    - (beta / 3) * ((1 + 6 * t + 9 * t.^2) / (1 + 4 * t).^(3/2) - 1) ...
    - (delta / 15) ...
        * ((1 + 10 * t + 30 * t.^2 + 25 * t.^3) ./ (1 + 4 * t).^(3/2) ...
            - (1 + 4 * t));
sddotfun = @(t) sddot0fun(t) ...
    - 2 * beta * (1 + 3 * t) ./ (1 + 4 * t).^(5/2) ...
    - (delta / 15) ...
        * ((4 + 40 * t + 135 * t.^2 + 150 * t.^3) ./ (1 + 4 * t).^(5/2) - 4);

% Solutions for s, sdot and sddot as arrays over the given time
s = sfun(tvals);
sdot = sdotfun(tvals);
sddot = sddotfun(tvals);

end
 
