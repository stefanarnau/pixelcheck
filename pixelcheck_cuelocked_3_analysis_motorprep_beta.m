clear all

% Paths 
PATH_EEGLAB = '/home/plkn/eeglab2025.0.0/';
PATH_TF_DATA = '/mnt/data_dump/pixelcheck/3_tf_data_cuelocked/';
%PATH_TF_DATA = '/mnt/data_dump/pixelcheck/3_tf_data_cuelocked_1stout/';
PATH_TF_RESULTS = '/mnt/data_dump/pixelcheck/4_tf_results_cuelocked/';

% Init EEGlab
addpath(PATH_EEGLAB);
eeglab;

% Get list of files
file_list = dir(fullfile(PATH_TF_DATA, '*ersps.mat'));
  
% Exclude VP4 VP25
to_exclude = [4, 25];
idx_exclude = [];
for s = 1 : length(file_list)

    % Get id as string
    id_string = regexp(file_list(s).name, '\d+', 'match');
    
    if ismember(str2double(id_string{1}), to_exclude)
        idx_exclude(end + 1) = s;
    end
end
file_list(idx_exclude) = [];

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

% Set label for analysis
analysis_label = 'motorprep_beta';

% Posterior electrode patch
idx_channel = [5, 6, 43, 44]; % motor

idx_time = tf_times >= 1200 & tf_times <= 1500;
idx_freq = tf_freqs >= 16 & tf_freqs <= 30;

% Average across electrodes in patch
ersps_patch = squeeze(mean(all_ersps(:, :, idx_channel, :, :), 3));

% Calculate values
ersps_vals = squeeze(mean(ersps_patch(:, :, idx_freq, idx_time), [3, 4]));

% Calculate freqtraces in patch
freq_patch = squeeze(mean(ersps_patch(:, :, idx_freq, :), 3));

% Get topovals for conditions
topovals = squeeze(mean(all_ersps(:, :, :, idx_freq, idx_time), [1, 4, 5]));

% Plot freqtraces
figure()
pd = squeeze(mean(freq_patch, 1));
plot(tf_times, pd)
legend({'neulo', 'neuhi', 'slflo', 'slfhi', 'othlo', 'othhi'})

tmp = tf_times(idx_time);
xRegion = [tmp(1), tmp(end), tmp(end), tmp(1)];
yRegion = [min(pd(:)), min(pd(:)), max(pd(:)), max(pd(:))];
h = patch(xRegion, yRegion, 'black');
set(h, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% Plot topo
plotlims = [-max(abs(topovals(:))), max(abs(topovals(:)))];
figure()
subplot(3, 2, 1)
topoplot(topovals(1, :), chanlocs, 'maplimits', plotlims, 'electrodes', 'on', 'style', 'both', 'emarker2', {idx_channel, 'o', 'r', 3, 2});
title(['neu lo'])
subplot(3, 2, 2)
topoplot(topovals(2, :), chanlocs, 'maplimits', plotlims, 'electrodes', 'on', 'style', 'both', 'emarker2', {idx_channel, 'o', 'r', 3, 2});
title(['neu hi'])

subplot(3, 2, 3)
topoplot(topovals(3, :), chanlocs, 'maplimits', plotlims, 'electrodes', 'on', 'style', 'both', 'emarker2', {idx_channel, 'o', 'r', 3, 2});
title(['slf lo'])
subplot(3, 2, 4)
topoplot(topovals(4, :), chanlocs, 'maplimits', plotlims, 'electrodes', 'on', 'style', 'both', 'emarker2', {idx_channel, 'o', 'r', 3, 2});
title(['slf hi'])

subplot(3, 2, 5)
topoplot(topovals(5, :), chanlocs, 'maplimits', plotlims, 'electrodes', 'on', 'style', 'both', 'emarker2', {idx_channel, 'o', 'r', 3, 2});
title(['oth lo'])
subplot(3, 2, 6)
topoplot(topovals(6, :), chanlocs, 'maplimits', plotlims, 'electrodes', 'on', 'style', 'both', 'emarker2', {idx_channel, 'o', 'r', 3, 2});
title(['oth hi'])


% Plot tf for conditions
figure()

subplot(3, 2, 1)
pd = squeeze(mean(squeeze(ersps_patch(:, 1, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'linecolor','none')
colormap('jet')
set(gca,'clim', [-2, 2], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
colorbar;
title(['neu lo'])

subplot(3, 2, 2)
pd = squeeze(mean(squeeze(ersps_patch(:, 2, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'linecolor','none')
colormap('jet')
set(gca,'clim', [-2, 2], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
colorbar;
title(['neu hi'])

subplot(3, 2, 3)
pd = squeeze(mean(squeeze(ersps_patch(:, 3, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'linecolor','none')
colormap('jet')
set(gca,'clim', [-2, 2], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
colorbar;
title(['slf lo'])

subplot(3, 2, 4)
pd = squeeze(mean(squeeze(ersps_patch(:, 4, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'linecolor','none')
colormap('jet')
set(gca,'clim', [-2, 2], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
colorbar;
title(['slf hi'])

subplot(3, 2, 5)
pd = squeeze(mean(squeeze(ersps_patch(:, 5, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'linecolor','none')
colormap('jet')
set(gca,'clim', [-2, 2], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
colorbar;
title(['oth lo'])

subplot(3, 2, 6)
pd = squeeze(mean(squeeze(ersps_patch(:, 6, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'linecolor','none')
colormap('jet')
set(gca,'clim', [-2, 2], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
colorbar;
title(['oth hi'])


% Create output matrix and save (long format)
out_long = [];
counter = 0;
for s = 1 : size(ersps_vals, 1)
    for cond = 1 : size(ersps_vals, 2)
        counter = counter + 1;
        if cond <= 2
            fb = 1;
        elseif cond <= 4
            fb = 2;
        else
            fb = 3;
        end
        if mod(cond, 2) == 1
            rew = 0;
        else
            rew = 1;
        end
        out_long(counter, :) = [s, fb, rew, ersps_vals(s, cond)];
 
    end
end
fn = [PATH_TF_RESULTS, analysis_label, '.csv'];
writematrix(out_long, fn);