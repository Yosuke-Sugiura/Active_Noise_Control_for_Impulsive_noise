function [w, params] = fxlognlms(w, e, r_buf, ~, mu, params)
% FxlogNLMS algorithm (Logarithmic-transformed FxNLMS)
%
% Paper:
% L. Wu, H. He, and X. Qiu, “An active impulsive noise control algorithm with logarithmic transformation,” 
% IEEE Transactions on Audio, Speech, and Language Processing, vol. 19, no. 4, pp. 1041–1044, 2011.
%
% Params:
%   G_l : Gain factor for logarithmic transformation

    sgn = sign(e);
    e_log = log(params.G_l * abs(e));
    if e_log < 0, e_log = 0; end
    e_c = sgn * e_log / (params.G_l * abs(e));
    N = length(r_buf);
    mu_ = mu / (N * mean(r_buf.^2) + 1e-3);
    w = w - mu_ * e_c * r_buf;
end
