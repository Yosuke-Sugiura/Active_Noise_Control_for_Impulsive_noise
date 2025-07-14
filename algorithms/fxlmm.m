function [w, params] = fxlmm(w, e, r, ~, mu, params)
% FxLMM (Hampel 3-part M-estimator, conventional version)
% Hampel break-points: delta1 = 3*delta0, delta2 = 4*delta0
%
% Paper:
% T. Thanigai, S. M. Kuo, and R. Yenduri, “Nonlinear active noise control for infant incubators in neo-natal intensive care units,” 
% in Proc. IEEE Int. Conf. Acoustics, Speech, and Signal Processing (ICASSP), vol. 1, 2007, pp. 109–112.
%
% Params:
%   params.gzai : Central threshold (delta0), controls the width of the linear region
%


    delta0 = params.gzai;          % central region width
    delta1 = 3 * delta0;           % first breakpoint
    delta2 = 4 * delta0;           % second breakpoint

    abs_e = abs(e);
    if abs_e < delta0
        psi = e;
    elseif abs_e < delta1
        psi = delta0 * sign(e);
    elseif abs_e < delta2
        psi = delta0 * sign(e) * (delta2 - abs_e) / (delta2 - delta1);
    else
        psi = 0;
    end

    w = w - mu * psi * r;
end
