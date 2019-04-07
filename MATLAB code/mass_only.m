%% Droplet impact on a cantilever: mass_only.m
% 
% This script deals with the distinguished limit where the mass term of the
% ODE for the cantilever dominates over the damping and spring term. In
% this case we have an analytic solution for s_0(t) for all time. We drop
% the 0 subscripts.

%%
% Solution for s_0(t) and other quantities
s = @(t) 0.5 + t - 0.5 * sqrt(1 + 4 * t);
sdot = @(t) 1 - 1 / sqrt(1 + 4 * t);
d = @(t) 2 * sqrt(t - s(t));
