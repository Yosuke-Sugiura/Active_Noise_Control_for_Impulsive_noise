function [w, params] = fair(w, e, r_buf, ~, mu, params)
% Fair algorithm (M-estimator based)
%
% Paper:
% L. Wu and X. Qiu, “An M-estimator based algorithm for active impulse-like noise control,” 
% Applied Acoustics, vol. 74, pp. 407–412, 2013.
%
% Params:
%   (none required explicitly; 'params' is included for consistency or future use)

    A = mean(abs(e));
    c = A;
    e_c = e / (1 + abs(e) / c);
    w = w - mu * e_c * r_buf;
end