clear all

% Paths 
PATH_EEGLAB = '/home/plkn/eeglab2025.0.0/';
PATH_TF_DATA = '/mnt/data_dump/pixelcheck/3_tf_data_cuelocked/';
PATH_TF_RESULTS = '/mnt/data_dump/pixelcheck/4_tf_results_cuelocked/';

% Init EEGlab
addpath(PATH_EEGLAB);
eeglab;

% Get list of files
file_list = dir(fullfile(PATH_TF_DATA, '*ersps.mat'));
  
% Get number of subjects
n_subjects = length(file_list);

% Load metadata
load([PATH_TF_DATA, 'chanlocs.mat']);
load([PATH_TF_DATA, 'tf_freqs.mat']);
load([PATH_TF_DATA, 'tf_times.mat']);

% Prune to cti
time_idx = tf_times >= -500 & tf_times <= 1800;
tf_times = tf_times(time_idx);

% Data matrix
all_ersps = zeros(n_subjects, 6, length(chanlocs), length(tf_freqs), length(tf_times));

% Loop subjects and load data
for s = 1 : n_subjects

    % Get id as string
    id_string = regexp(file_list(s).name, '\d+', 'match');
    id_string = id_string{1};

    % Load data
    load([PATH_TF_DATA, id_string, '_ersps.mat']); % cond x chan x freq x time

    % Prune to cti
    ersps = ersps(:, :, :, time_idx);

    % Collect
    all_ersps(s, :, :, :, :) = ersps;
   
end

% Posterior electrode patch
idx_posterior = [17, 18, 64];
ersps_posterior = squeeze(mean(all_ersps(:, :, idx_posterior, :, :), 3));


figure()

subplot(3, 2, 1)
pd = squeeze(mean(squeeze(ersps_posterior(:, 1, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'linecolor','none')
colormap('jet')
set(gca,'clim', [-2, 2], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
colorbar;
title(['neu lo'])

subplot(3, 2, 2)
pd = squeeze(mean(squeeze(ersps_posterior(:, 2, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'linecolor','none')
colormap('jet')
set(gca,'clim', [-2, 2], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
colorbar;
title(['neu hi'])

subplot(3, 2, 3)
pd = squeeze(mean(squeeze(ersps_posterior(:, 3, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'linecolor','none')
colormap('jet')
set(gca,'clim', [-2, 2], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
colorbar;
title(['slf lo'])

subplot(3, 2, 4)
pd = squeeze(mean(squeeze(ersps_posterior(:, 4, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'linecolor','none')
colormap('jet')
set(gca,'clim', [-2, 2], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
colorbar;
title(['slf hi'])

subplot(3, 2, 5)
pd = squeeze(mean(squeeze(ersps_posterior(:, 5, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'linecolor','none')
colormap('jet')
set(gca,'clim', [-2, 2], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
colorbar;
title(['oth lo'])

subplot(3, 2, 6)
pd = squeeze(mean(squeeze(ersps_posterior(:, 6, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'linecolor','none')
colormap('jet')
set(gca,'clim', [-2, 2], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
colorbar;
title(['oth hi'])



