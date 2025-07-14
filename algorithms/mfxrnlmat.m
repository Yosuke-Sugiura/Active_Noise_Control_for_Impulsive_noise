function [w, params] = mfxrnlmat(w, e, r_buf, ~, mu, params)
% MFxRNLMAT algorithm (Robust Nonlinear Modified Adaptive Technique)
%
% Paper:
% A. Mirza, F. Afzal, A. Zeb, A. Wakeel, W. S. Qureshi, and A. Akgul, 
% “New FxLMAT-Based Algorithms for Active Control of Impulsive Noise,” 
% IEEE Access, vol. 11, pp. 81279–81288, 2023, doi: 10.1109/ACCESS.2023.3293647.
%
% Params:
%   lmd_m  : Forgetting factor for error energy smoothing (Ke)
%   Ke     : Smoothed squared error energy (updated internally)
%   dlt2   : Regularization constant
%   Bias   : Offset bias to stabilize denominator
%   beta2  : Gain factor for nonlinear normalization term
    
    params.Ke = params.lmd_m * params.Ke + (1 - params.lmd_m) * abs(e)^2;
    
    mu_m = mu / (sum(r_buf.^2) + params.Ke + params.dlt2);
    mu_  = mu_m / (1 + params.beta2 * abs(e)^3);
    
    w = w - mu_ * abs(e)^2 * sign(e) * r_buf;
end
