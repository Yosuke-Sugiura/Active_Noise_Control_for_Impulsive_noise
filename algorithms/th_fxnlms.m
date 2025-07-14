function [w, params] = th_fxnlms(w, e, r_buf, ~, mu, params)
    % Thresholded FxNLMS
    x_c = min(max(e, params.c1), params.c2);
    N = length(r_buf);
    mu_ = mu / (N * mean(r_buf.^2) + 1e-3);
    w = w - mu_ * x_c * r_buf;
end