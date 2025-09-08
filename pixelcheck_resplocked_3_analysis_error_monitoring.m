clear all

% Paths 
PATH_EEGLAB = '/home/plkn/eeglab2025.0.0/';
PATH_FIELDTRIP = '/home/plkn/fieldtrip-master/';
PATH_TF_DATA = '/mnt/data_dump/pixelcheck/3_tf_data_resplocked/';
PATH_TF_RESULTS = '/mnt/data_dump/pixelcheck/4_tf_results/';

% Init EEGlab
addpath(PATH_EEGLAB);
eeglab;

% Init fieldtrip
addpath(PATH_FIELDTRIP);
ft_defaults;

% Get list of files
file_list = dir(fullfile(PATH_TF_DATA, '*ersps.mat'));
    
% Remember condition labels
condizion_labels = {'slf err lo',...
                    'slf err hi',...
                    'slf flip lo',...
                    'slf flip hi',...
                    'slf corr lo',... 
                    'slf corr hi',...
                    'oth err lo',...
                    'oth err hi',...
                    'oth flip lo',...
                    'oth flip hi',...
                    'oth corr lo',...
                    'oth corr hi',...
                    'neu err lo',...
                    'neu err hi',...
                    'neu corr lo',...
                    'neu corr hi'};
                
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

% Data matrix
error_diffs = zeros(n_subjects, 6, length(chanlocs), length(tf_freqs), length(tf_times));

% Loop subjects
for s = 1 : n_subjects

    % Get participant id
    id_string = regexp(file_list(s).name, '\d+', 'match'); 

    % Load
    load([PATH_TF_DATA, id_string{1}, '_ersps.mat']); % cond x chan x freq x time

    error_diffs(s, 1, :, :, :) = squeeze(ersps(13, :, :, :));
    error_diffs(s, 2, :, :, :) = squeeze(ersps(14, :, :, :));
    error_diffs(s, 3, :, :, :) = squeeze(ersps(1, :, :, :));
    error_diffs(s, 4, :, :, :) = squeeze(ersps(2, :, :, :));
    error_diffs(s, 5, :, :, :) = squeeze(ersps(7, :, :, :));
    error_diffs(s, 6, :, :, :) = squeeze(ersps(8, :, :, :));

end % end subject loop

% Set label for analysis
analysis_label = 'error_monitoring_frontal_theta';

% Posterior electrode patch
idx_channel = [65, 15, 19, 20, 16]; % frontal midline

idx_time = tf_times >= 50 & tf_times <= 250;
idx_freq = tf_freqs >= 4 & tf_freqs <= 7;

% Average across electrodes in patch
ersps_patch = squeeze(mean(error_diffs(:, :, idx_channel, :, :), 3));

% Calculate values
ersps_vals = squeeze(mean(ersps_patch(:, :, idx_freq, idx_time), [3, 4]));

% Calculate freqtraces in patch
freq_patch = squeeze(mean(ersps_patch(:, :, idx_freq, :), 3));

% Get topovals for conditions
topovals = squeeze(mean(error_diffs(:, :, :, idx_freq, idx_time), [1, 4, 5]));

% Plot freqtraces
figure()
pd = squeeze(mean(freq_patch, 1));
plot(tf_times, pd)
legend({'neu lo', 'neu hi', 'slf lo', 'slf hi','oth lo', 'oth hi'})
tmp = tf_times(idx_time);
xRegion = [tmp(1), tmp(end), tmp(end), tmp(1)];
yRegion = [min(pd(:)), min(pd(:)), max(pd(:)), max(pd(:))];
h = patch(xRegion, yRegion, 'black');
set(h, 'FaceAlpha', 0.1, 'EdgeColor', 'none');


% Plot topo
plotlims = [-max(abs(topovals(:))), max(abs(topovals(:)))];
figure()
subplot(3, 2, 1)
topoplot(topovals(1, :), chanlocs, 'maplimits', plotlims, 'electrodes', 'on', 'style', 'both', 'emarker2', {idx_channel, 'o', 'c', 3, 2});
title(['neu lo'])

subplot(3, 2, 2)
topoplot(topovals(2, :), chanlocs, 'maplimits', plotlims, 'electrodes', 'on', 'style', 'both', 'emarker2', {idx_channel, 'o', 'c', 3, 2});
title(['neu hi'])

subplot(3, 2, 3)
topoplot(topovals(3, :), chanlocs, 'maplimits', plotlims, 'electrodes', 'on', 'style', 'both', 'emarker2', {idx_channel, 'o', 'c', 3, 2});
title(['slf lo'])

subplot(3, 2, 4)
topoplot(topovals(4, :), chanlocs, 'maplimits', plotlims, 'electrodes', 'on', 'style', 'both', 'emarker2', {idx_channel, 'o', 'c', 3, 2});
title(['slf hi'])

subplot(3, 2, 5)
topoplot(topovals(5, :), chanlocs, 'maplimits', plotlims, 'electrodes', 'on', 'style', 'both', 'emarker2', {idx_channel, 'o', 'c', 3, 2});
title(['oth lo'])

subplot(3, 2, 6)
topoplot(topovals(6, :), chanlocs, 'maplimits', plotlims, 'electrodes', 'on', 'style', 'both', 'emarker2', {idx_channel, 'o', 'c', 3, 2});
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
    for cond = 1 : 6
        counter = counter + 1;
        if cond == 1
            val = ersps_vals(s, 1);
            fb = 1; % slf
            rew = 1; % lo
        elseif cond == 2
            val = ersps_vals(s, 2);
            fb = 1; % slf
            rew = 2; % hi
        elseif cond == 3
            val = ersps_vals(s, 3);
            fb = 2; % slf
            rew = 1; % lo
        elseif cond == 4
            val = ersps_vals(s, 4);
            fb = 2; % slf
            rew = 2; % hi
        elseif cond == 5
            val = ersps_vals(s, 5);
            fb = 3; % oth
            rew = 1; % lo
        elseif cond == 6
            val = ersps_vals(s, 6);
            fb = 3; % oth
            rew = 2; % hi
        end
        out_long(counter, :) = [s, fb, rew, val];
    end
end
fn = [PATH_TF_RESULTS, analysis_label, '.csv'];
writematrix(out_long, fn);