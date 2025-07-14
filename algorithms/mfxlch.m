function [w, params] = mfxlch(w, e, r_buf, ~, mu, params)
% MFxLCH algorithm (Modified Filtered-x Least-squares Correntropy-based Hybrid)
%
% Paper:
% A. Mirza, A. Zeb, M. Y. Umair, D. Ilyas, and S. A. Sheikh, “Less complex solutions for active noise control of impulsive noise,” 
% Analog Integrated Circuits and Signal Processing, vol. 102, no. 3, pp. 507–521, Mar. 2020.
%
% Params:
%   rho    : Threshold parameter for nonlinear error weighting
%   Ke     : Energy smoothing term (updated internally)
%   lmd_m  : Forgetting factor for Ke update
%   dlt2   : Regularization term to prevent division by zero

    rho = params.rho;
    if abs(e) > 1/rho
        e_c = sign(e);
    else
        e_c = -e * abs(e) * rho^2 + 2 * rho * e;
    end

    params.Ke = params.lmd_m * params.Ke + (1 - params.lmd_m) * e^2;
    mu_ = mu / (length(r_buf) * mean(r_buf.^2) + params.Ke + params.dlt2);
    w = w - mu_ * e_c * r_buf;
end