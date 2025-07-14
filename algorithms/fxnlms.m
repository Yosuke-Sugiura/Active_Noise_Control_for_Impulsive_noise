function [w, params] = fxnlms(w, e, r_buf, ~, mu, params)
% FxNLMS algorithm
    mu_ = mu / (sum(r_buf.^2) + 1e-3);
    w = w - mu_ * e * r_buf;
end