function [w, params] = fxlognlms_plus(w, e, r_buf, ~, mu, params)
% NSS-FxlogLMS+ algorithm (Nonlinear Step-size Normalized FxlogLMS Plus)
%
% Paper:
% A. Haneda, Y. Sugiura, and T. Shimamura, “FxlogLMS+: Modified FxlogLMS Algorithm for Active Impulsive Noise Control,” 
% in Proc. ICGEC 2024: Genetic and Evolutionary Computing, Lecture Notes in Electrical Engineering, vol. 1322, 
% Springer, Singapore, 2025, pp. 342-351, doi: 10.1007/978-981-96-1535-3_34.
%
% Params:
%   G_ml : Gain factor for error shaping

    sgn = sign(e);
    abs_e = params.G_ml * abs(e) + 1;
    e_c = sgn * log(abs_e) / abs_e;
    N = length(r_buf);
    mu_ = mu / (N * mean(r_buf.^2) + 1e-3);
    w = w - mu_ * e_c * r_buf;
end