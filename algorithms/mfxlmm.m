function [w, params] = mfxlmm(w, e, r_buf, ~, mu, params)
% M‑FxLMM  (Modified Filtered‑x Least‑Mean‑M‑estimate)
%
% Paper:
% G. Sun, M. T. Li, and T. C. Lim, “Enhanced filtered‑x least mean M‑estimate algorithm for active impulsive noise control,” *Applied Acoustics*, vol. 90, pp. 31–41, Oct. 2015, doi: 10.1016/j.apacoust.2015.05.007.
%
% params:
% - gzai_l   : Hampel δ₀
% - sigma_l  : initial σ_e²
% - lmd_l    : forgetting factor for σ_e² (0.9–0.99)
% - eps_l    : small ε > 0 to avoid zero‑division


    % ----- Hampel score -----
    d0 = params.gzai_l;          % δ₀
    d1 = 3 * d0;               % δ₁
    d2 = 4 * d0;               % δ₂

    abs_e = abs(e);
    if abs_e < d0
        psi = e;
    elseif abs_e < d1
        psi = d0 * sign(e);
    elseif abs_e < d2
        psi = d0 * sign(e) * (d2 - abs_e) / (d2 - d1);
    else
        psi = 0;
    end

    % ----- update σ_e² (error energy) -----
    params.sigma_l = params.lmd_l * params.sigma_l + ...
                     (1 - params.lmd_l) * abs(e)^2;

    % ----- adaptive step -----
    power_r = mean(r_buf.^2);
    mu_adapt = mu / (power_r + params.sigma_l + params.eps_l);

    % ----- weight update -----
    w = w - mu_adapt * psi * r_buf;
end
