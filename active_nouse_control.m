% Active Noise Control - Feedforward FxLMS Variants
% Author: Yosuke Sugiura, Aoi Hada (Saitama University, Japan)
% -------------------------------------------------
% Ver. 1    Date: 2024-04-18

clear; close all;
addpath('algorithms'); % Load algorithm functions

%% Parameters
fs = 16000;                 % Sampling frequency
order_control = 512;        % Order of control filter W(z)
order_secondary = 128;      % Order of secondary path model C_h(z)
frame_step = 64;            % Frame step size for ANR
num_trials = 1;             % Number of trials

%% Load noise signal
% Select a noise file from the following options:
%  "sounds/impulsive_noise_alpha1.65.dat" : alpha = 1.65  (30s, 16kHz)
%  "sounds/impulsive_noise_alpha1.45.dat" : alpha = 1.45  (30s, 16kHz)
%  "sounds/impulsive_noise_alpha1.85.dat" : alpha = 1.85  (30s, 16kHz)
%  "sounds/impulsive_noise_alpha2.00.dat" : alpha = 2.00  (30s, 16kHz)
%  "sounds/impulsive_noise_alpha1.65_1.45_1.85.dat" : alpha = 1.65 (0~10s),1.45 (10s~20s), 1.85 (20s~30s, 16kHz)
file = "sounds/impulsive_noise_alpha1.45.dat";
noise = load(file);
signal_length = length(noise);

%% Algorithm setup
% alg_names : Please choose algorithms to run from...
%   "FxNLMS", "FxgsnLMS", "FxLMM", "MFxLMM", "Th_FxNLMS", "MFxLCH", "MFxRNLMAT", "Fair", "FxlogNLMS", "NS_FxlogLMS", "FxlogNLMS_plus", "NSS_FxlogLMS_plus"
alg_names = ["FxNLMS", "FxgsnLMS", "FxLMM", "MFxLMM", "Th_FxNLMS", "MFxLCH", "MFxRNLMAT", "Fair", "FxlogNLMS", "NS_FxlogLMS", "FxlogNLMS_plus", "NSS_FxlogLMS_plus"];
plot_colors = {'#000000', '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#bcbd22', '#17becf', '#9edae5', '#c49c94', '#f7b6d2', '#c5b0d5', '#aec7e8'};

mu_values = dictionary(...
    "FxNLMS", 0.25, "FxgsnLMS", 0.5, "FxLMM", 0.001, "MFxLMM", 0.001, "MFxLCH", 0.03, "MFxRNLMAT", 15.0, "Th_FxNLMS", 0.3, ...
    "Fair", 0.001, "FxlogNLMS", 0.5, "NS_FxlogLMS", 0.1, "FxlogNLMS_plus", 0.1, "NSS_FxlogLMS_plus", 0.1);

params_template = struct(...
    'c1', -200, 'c2', 200, 'sigma_e', 1, 'G_l', 1e4, 'lmd_e', 0.99, 'G_ml', 5, ...
    'rho', 10, 'lmd_m', 0.9999, 'dlt2', 0.01, 'beta2', 5000, 'Bias', 0.001, 'gzai', 5e-2, 'gzai_l', 5e-2, ...
    'cx1', -50, 'cx2', 50, 'ce1', -50, 'ce2', 50, 'rho1eps1', 0.01, 'rho2', 0.01, 'eps2', 0.01, ...
    'sigma_l',1, 'lmd_l',0.97, 'eps_l',1e-3, 'NS_t', 1.0, 'E_e', 0, 'Ke', 0);

%% Load impulse responses
imp_primary = load('impulse_response/primary_ir.dat');
imp_secondary = load('impulse_response/secondary_ir.dat');
imp_primary = [imp_primary; zeros(order_control,1)];
imp_secondary = [imp_secondary; zeros(order_secondary,1)];
imp_primary = imp_primary(1:order_control)';
imp_secondary = imp_secondary(1:order_secondary)';

%% Visualize impulse response characteristics (both paths in one figure)
freq_axis = (0:511) * fs / 512;
fft_primary = fft(imp_primary, 512);
fft_secondary = fft(imp_secondary, 512);

figure('Name', 'Impulse Response Characteristics');
subplot(2,1,1);
plot(freq_axis, 20*log10(abs(fft_primary)+1e-2), 'b', 'DisplayName', 'Primary'); hold on;
plot(freq_axis, 20*log10(abs(fft_secondary)+1e-2), 'r', 'DisplayName', 'Secondary');
xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]'); title('Magnitude Response'); xlim([0 fs/2]); legend;

subplot(2,1,2);
plot(freq_axis, unwrap(angle(fft_primary)), 'b', 'DisplayName', 'Primary'); hold on;
plot(freq_axis, unwrap(angle(fft_secondary)), 'r', 'DisplayName', 'Secondary');
xlabel('Frequency [Hz]'); ylabel('Phase [rad]'); title('Phase Response'); xlim([0 fs/2]); legend;

%% Initialize buffers
ref_signal = zeros(signal_length,1);
error_signal = zeros(signal_length,length(alg_names));
avg_ANR = zeros(floor(signal_length/frame_step),length(alg_names));

%% Run simulation
fprintf('Process start : %s\n', file)
for alg_idx = 1:length(alg_names)
    alg_name = alg_names{alg_idx};
    mu = mu_values(alg_name);
    params = params_template;
    update_func = str2func(lower(alg_name));

    for trial = 1:num_trials
        input_noise = noise(:,trial);

        w = zeros(1, order_control);
        ch = imp_secondary(1:order_secondary);
        buf_x = zeros(max(order_control, order_secondary),1);
        buf_y = zeros(max(length(imp_secondary),order_secondary),1);
        buf_r = zeros(1, order_control);
        anr_temp = zeros(floor(signal_length/frame_step),1);
        e_m = 0; d_m = 0; n_frame = 1;

        for n = 1:signal_length
            x = input_noise(n);
            buf_x = [x; buf_x(1:end-1)];

            y_primary = imp_primary * buf_x(1:order_control);
            y_control = w * buf_x(1:order_control);
            buf_y = [y_control; buf_y(1:end-1)];
            y_secondary = imp_secondary * buf_y(1:length(imp_secondary));

            r = ch * buf_x(1:order_secondary);
            buf_r = [r, buf_r(1:end-1)];

            e = y_primary + y_secondary;
            [w, params] = update_func(w, e, buf_r, buf_x, mu, params);

            if rem(n, frame_step) == 0
                anr_temp(n_frame) = 20 * log10(e_m / d_m + 1e-8);
                n_frame = n_frame + 1;
            end

            ref_signal(n) = y_primary;
            error_signal(n, alg_idx) = e;
            e_m = 0.999 * e_m + 0.001 * abs(e);
            d_m = 0.999 * d_m + 0.001 * abs(y_primary);
        end

        avg_ANR(:,alg_idx) = avg_ANR(:,alg_idx) + anr_temp / num_trials;
        fprintf(' [%s]   Average ANR = %.2f dB\n', alg_name, mean(avg_ANR(:,alg_idx)));
    end
end

%% Plot ANR
figure('Name', 'ANR'); hold on;
time_axis = (0:floor(signal_length/frame_step)-1) * frame_step / fs;
for i = 1:length(alg_names)
    plot(time_axis, avg_ANR(:,i), 'Color', plot_colors{i}, 'DisplayName', strrep(alg_names{i}, '_', '-'));
end
xlabel('Time [s]'); ylabel('ANR (dB)'); legend;

%% Save results
save('logs/input.dat', 'ref_signal', '-ascii');
for i = 1:length(alg_names)
    out_i = error_signal(:,i);
    save(['logs/output_', alg_names{i}, '.dat'], 'out_i', '-ascii');
end