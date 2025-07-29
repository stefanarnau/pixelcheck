clear all;

% File path
PATH_EEGLAB = '/home/plkn/eeglab2025.0.0/';
PATH_IN = '/mnt/data_dump/pixelcheck_resplock_ged/ged_calculated_flipped_theta/';

% Init EEGlab
addpath(PATH_EEGLAB);
eeglab;

% Get list of files
file_list = dir(fullfile(PATH_IN, '*.mat'));

% Load first dataset
fn = file_list(1).name;
load([PATH_IN, fn])

% Wavelet Parameters
fs = 200;
n_frq = 30;
frqrange = [2, 30];
tfres_range = [600, 240];
wtime = -2 : 1 / fs : 2;

hz = linspace(0, fs, length(wtime));

tf_freqs = logspace(log10(frqrange(1)), log10(frqrange(2)), n_frq);
fwhmTs = logspace(log10(tfres_range(1)), log10(tfres_range(2)), n_frq);

% Create wavelets
cmw = zeros(length(tf_freqs), length(wtime));
for frq = 1:length(tf_freqs)
    cmw(frq,:) = exp(2*1i*pi*tf_freqs(frq).*wtime) .* exp((-4*log(2)*wtime.^2) ./ (fwhmTs(frq)/1000)^2);
    cmw(frq,:) = cmw(frq,:) ./ max(cmw(frq,:));
end

% Crop times
prune_times = [-500, 1500];
tf_times = ged_data.times(dsearchn(ged_data.times', prune_times(1)):dsearchn(ged_data.times', prune_times(2)));

% Result matrices
ersps = zeros(length(file_list), 12, length(tf_freqs), length(tf_times));
itpcs = zeros(length(file_list), 12, length(tf_freqs), length(tf_times));

% Iterate files
for file_number = 1 : length(file_list)

    % Get filename
    fn = file_list(file_number).name;

    % Talk
    fprintf('\ntf-decomposing %s\n', fn);

    % Load
    load([PATH_IN, fn])

    % Code trialnumber in block
    ged_data.trialinfo.trial_in_block = mod(ged_data.trialinfo.trial_nr, 120);
    ged_data.trialinfo.trial_in_block(ged_data.trialinfo.trial_in_block == 0) = 120;

    % Set trial exclusion criteria
    idx_keep = ged_data.trialinfo.trial_in_block > 40 & ged_data.trialinfo.block_nr >=2;

    % Exclude trials
    signal = ged_data.ged_time_series(:, idx_keep);
    trialinfo = ged_data.trialinfo(idx_keep, :);

    % Perform fft on signal
    convlen = size(signal, 1) * size(signal, 2) + size(cmw, 2) - 1;
    tmp = fft(reshape(signal, 1, []), convlen);

    powcube = NaN(length(tf_freqs), size(signal, 1), size(signal, 2));
    phacube = NaN(length(tf_freqs), size(signal, 1), size(signal, 2));

    % Convolute
    cmwX = zeros(length(tf_freqs), convlen);
    for f = 1 : length(tf_freqs)
        cmwX(f, :) = fft(cmw(f, :), convlen);
        cmwX(f, :) = cmwX(f, :) ./ max(cmwX(f, :));
    end

    for f = 1 : length(tf_freqs)
        as = ifft(cmwX(f, :) .* tmp);
        as = as(((size(cmw, 2) - 1) / 2) + 1 : end - ((size(cmw, 2) - 1) / 2));
        as = reshape(as, size(signal, 1), size(signal, 2));
        powcube(f, :, :) = abs(as) .^ 2;
        phacube(f, :, :) = exp(1i * angle(as));
    end

    % Crop time
    idx_t1 = dsearchn(ged_data.times', prune_times(1));
    idx_t2 = dsearchn(ged_data.times', prune_times(2));
    powcube = powcube(:, idx_t1 : idx_t2, :);
    phacube = phacube(:, idx_t1 : idx_t2, :);

    % Baseline [-500:-200]
    [~, bl1] = min(abs(tf_times - (-500)));
    [~, bl2] = min(abs(tf_times - (-200)));
    blvals = squeeze(mean(mean(powcube(:, bl1 : bl2, :), 3), 2));

    % Loop conditions
    for cond = 1 : 12
        switch cond
            case 1 % self incorrect low
                idx = trialinfo.ma_condition==2 & trialinfo.fb_correct==0 & trialinfo.fb_flipped==0 & trialinfo.reward_condition==0;
            case 2 % self incorrect high
                idx = trialinfo.ma_condition==2 & trialinfo.fb_correct==0 & trialinfo.fb_flipped==0 & trialinfo.reward_condition==1;
            case 3
                idx = trialinfo.ma_condition==2 & trialinfo.fb_correct==0 & trialinfo.fb_flipped==1 & trialinfo.reward_condition==0;
            case 4
                idx = trialinfo.ma_condition==2 & trialinfo.fb_correct==0 & trialinfo.fb_flipped==1 & trialinfo.reward_condition==1;
            case 5
                idx = trialinfo.ma_condition==2 & trialinfo.fb_correct==1 & trialinfo.reward_condition==0;
            case 6
                idx = trialinfo.ma_condition==2 & trialinfo.fb_correct==1 & trialinfo.reward_condition==1;
            case 7
                idx = trialinfo.ma_condition==3 & trialinfo.fb_correct==0 & trialinfo.fb_flipped==0 & trialinfo.reward_condition==0;
            case 8
                idx = trialinfo.ma_condition==3 & trialinfo.fb_correct==0 & trialinfo.fb_flipped==0 & trialinfo.reward_condition==1;
            case 9
                idx = trialinfo.ma_condition==3 & trialinfo.fb_correct==0 & trialinfo.fb_flipped==1 & trialinfo.reward_condition==0;
            case 10
                idx = trialinfo.ma_condition==3 & trialinfo.fb_correct==0 & trialinfo.fb_flipped==1 & trialinfo.reward_condition==1;
            case 11
                idx = trialinfo.ma_condition==3 & trialinfo.fb_correct==1 & trialinfo.reward_condition==0;
            case 12
                idx = trialinfo.ma_condition==3 & trialinfo.fb_correct==1 & trialinfo.reward_condition==1;
        end

        % Save to matrix
        ersps(file_number, cond, :, :) = 10*log10(bsxfun(@rdivide, mean(powcube(:, :, idx), 3), blvals));
        itpcs(file_number, cond, :, :) = abs(mean(phacube(:, :, idx), 3));

    end

end



figure()
cond_order = [5,  1, 3,...
              6,  2, 4,...
              11, 7, 9,...
              12, 8, 10];

for cond = 1 : 12

    subplot(4, 3, cond_order(cond))
    pd = squeeze(mean(squeeze(ersps(:, cond, :, :)), 1));
    contourf(tf_times, tf_freqs, pd, 50, 'LineStyle', 'none')
    colormap('jet')
    set(gca, 'clim', [-3, 3], 'xlim', [-500, 1000], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
end



figure()
subplot(2, 2, 1)
pd = squeeze(mean(squeeze(ersps(:, 3, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'LineStyle', 'none')
colormap('jet')
set(gca, 'clim', [-3, 3], 'xlim', [-500, 1000], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
title('self lo')

subplot(2, 2, 2)
pd = squeeze(mean(squeeze(ersps(:, 4, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'LineStyle', 'none')
colormap('jet')
set(gca, 'clim', [-3, 3], 'xlim', [-500, 1000], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
title('self hi')

subplot(2, 2, 3)
pd = squeeze(mean(squeeze(ersps(:, 9, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'LineStyle', 'none')
colormap('jet')
set(gca, 'clim', [-3, 3], 'xlim', [-500, 1000], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
title('other lo')

subplot(2, 2, 4)
pd = squeeze(mean(squeeze(ersps(:, 10, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'LineStyle', 'none')
colormap('jet')
set(gca, 'clim', [-3, 3], 'xlim', [-500, 1000], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
title('other hi')


theta_idx = tf_freqs >= 8 & tf_freqs <= 12;

ersp_theta = squeeze(mean(ersps(:, :, theta_idx, :), 3));


theta_ave = squeeze(mean(ersp_theta, 1));

figure()
subplot(1, 2, 1)
plot(tf_times, theta_ave(3, :), "-k")
hold on
plot(tf_times, theta_ave(4, :), ":k")
plot(tf_times, theta_ave(9, :), "-r")
plot(tf_times, theta_ave(10, :), ":r")
ylim([-2, 3])
title('flip')

subplot(1, 2, 2)
plot(tf_times, theta_ave(1, :), "-k")
hold on
plot(tf_times, theta_ave(2, :), ":k")
plot(tf_times, theta_ave(7, :), "-r")
plot(tf_times, theta_ave(8, :), ":r")
ylim([-2, 3])
title('incorrect')


