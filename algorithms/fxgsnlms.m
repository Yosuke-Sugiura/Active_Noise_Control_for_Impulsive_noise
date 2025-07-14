function [w, params] = fxgsnlms(w, e, r_buf, x_buf, mu, params)
% FxgsnLMS algorithm
% Mechanical Systems and Signal Processing, vol. 56, pp. 320–333, May 2015, doi: 10.1016/j.ymssp.2014.10.002.
%
% Paper:
% Y. Zhou, Q. Zhang, and Y. Yin, “Active control of impulsive noise with symmetric α‑stable distribution based on an improved step‑size normalized adaptive algorithm,” Mechanical Systems and Signal Processing, vol. 56, pp. 320–333, May 2015, doi: 10.1016/j.ymssp.2014.10.002.
%
% params:
% - sigma_e: Estimated standard deviation of the error signal (used in the noise model)

    mu_ = mu / (sqrt(2*pi)*params.sigma_e) * exp(-(x_buf(1)^2) / (2*params.sigma_e^2));
    if mu_ >= 1/sum(x_buf.^2)
        mu_ = 1/sum(x_buf.^2);
    end
    w = w - mu_ * e * r_buf;
end