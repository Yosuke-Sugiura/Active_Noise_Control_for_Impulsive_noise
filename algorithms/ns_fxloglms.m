function [w, params] = ns_fxloglms(w, e, r_buf, ~, mu, params)
% NS-FxlogLMS algorithm (Nonlinear Selection FxlogLMS)
%
% Paper:
% M. Pawelczyk, W. Wierzchowski, L. Wu, and X. Qiu, “An extension to the filtered-x LMS algorithm with logarithmic transformation,” 
% in *Proc. IEEE Int. Symposium on Signal Processing and Information Technology (ISSPIT)*, Abu Dhabi, UAE, 2015, pp. 454–459, 
% doi: 10.1109/ISSPIT.2015.7394378.
%
% Params:
%   G_l  : Gain factor for logarithmic transformation
%   NS_t : Threshold for switching between linear and log transformation

    if abs(e) > params.NS_t
        sgn = sign(e);
        e_log = log(params.G_l * abs(e));
        e_c = sgn * e_log / (params.G_l * abs(e));
    else
        e_c = e;
    end

    N = length(r_buf);
    mu_ = mu / (N * mean(r_buf.^2) + 1e-3);
    w = w - mu_ * e_c * r_buf;
end