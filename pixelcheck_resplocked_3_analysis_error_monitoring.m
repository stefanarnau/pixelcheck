clear all

% Paths 
PATH_EEGLAB = '/home/plkn/eeglab2025.0.0/';
PATH_FIELDTRIP = '/home/plkn/fieldtrip-master/';
PATH_TF_DATA = '/mnt/data_dump/pixelcheck/3_tf_data_resplocked/';
%PATH_TF_DATA = '/mnt/data_dump/pixelcheck/3_tf_data_resplocked_1stout/';
PATH_TF_RESULTS = '/mnt/data_dump/pixelcheck/4_tf_results_resplocked/';

% Init EEGlab
addpath(PATH_EEGLAB);
eeglab;

% Init fieldtrip
addpath(PATH_FIELDTRIP);
ft_defaults;

% Get list of files
file_list = dir(fullfile(PATH_TF_DATA, '*ersps.mat'));
    
% Remember condition labels
cond_label = {'slf err lo',...
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
              'oth corr hi'};
                
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
all_ersps = zeros(n_subjects, 12, length(chanlocs), length(tf_freqs), length(tf_times));

% Loop subjects
for s = 1 : n_subjects

    % Get participant id
    id_string = regexp(file_list(s).name, '\d+', 'match'); 

    % Load
    load([PATH_TF_DATA, id_string{1}, '_ersps.mat']); % cond x chan x freq x time

    % Collect
    all_ersps(s, :, :, :, :) = ersps;

end % end subject loop

% Set label for analysis
analysis_label = 'error_monitoring_frontal_theta';

% Posterior electrode patch
idx_channel = [65, 15, 19, 20]; % frontal midline

idx_time = tf_times >= 100 & tf_times <= 200;
idx_freq = tf_freqs >= 4 & tf_freqs <= 7;

% Average across electrodes in patch
ersps_patch = squeeze(mean(all_ersps(:, :, idx_channel, :, :), 3));

% Calculate values
ersps_vals = squeeze(mean(ersps_patch(:, :, idx_freq, idx_time), [3, 4]));

% Calculate freqtraces in patch
freq_patch = squeeze(mean(ersps_patch(:, :, idx_freq, :), 3));

% Get topovals for conditions
topovals = squeeze(mean(all_ersps(:, :, :, idx_freq, idx_time), [1, 4, 5]));

% The condition labels
condizion_labels = {'slf err lo', 'slf err hi', 'slf flip lo', 'slf flip hi', 'slf corr lo', 'slf corr hi',...
                    'oth err lo', 'oth err hi', 'oth flip lo', 'oth flip hi', 'oth corr lo', 'oth corr hi'};

% Plot freqtraces
figure()
pd = squeeze(mean(freq_patch, 1));
plot(tf_times, pd)
legend(condizion_labels)
tmp = tf_times(idx_time);
xRegion = [tmp(1), tmp(end), tmp(end), tmp(1)];
yRegion = [min(pd(:)), min(pd(:)), max(pd(:)), max(pd(:))];
h = patch(xRegion, yRegion, 'black');
set(h, 'FaceAlpha', 0.1, 'EdgeColor', 'none');


% Plot topo
plotlims = [-max(abs(topovals(:))), max(abs(topovals(:)))];
figure()
subplot(3, 3, 1)
topoplot(topovals(1, :), chanlocs, 'maplimits', plotlims, 'electrodes', 'on', 'style', 'both', 'emarker2', {idx_channel, 'o', 'r', 3, 2});
title(['slf lo'])

subplot(3, 3, 2)
topoplot(topovals(2, :), chanlocs, 'maplimits', plotlims, 'electrodes', 'on', 'style', 'both', 'emarker2', {idx_channel, 'o', 'r', 3, 2});
title(['slf hi'])

subplot(3, 3, 4)
topoplot(topovals(7, :), chanlocs, 'maplimits', plotlims, 'electrodes', 'on', 'style', 'both', 'emarker2', {idx_channel, 'o', 'r', 3, 2});
title(['oth lo'])

subplot(3, 3, 5)
topoplot(topovals(8, :), chanlocs, 'maplimits', plotlims, 'electrodes', 'on', 'style', 'both', 'emarker2', {idx_channel, 'o', 'r', 3, 2});
title(['oth hi'])

subplot(3, 3, 3)
topoplot(topovals(2, :) - topovals(1, :), chanlocs, 'maplimits', plotlims, 'electrodes', 'on', 'style', 'both', 'emarker2', {idx_channel, 'o', 'r', 3, 2});
title(['slf hi-lo'])

subplot(3, 3, 6)
topoplot(topovals(8, :) - topovals(7, :), chanlocs, 'maplimits', plotlims, 'electrodes', 'on', 'style', 'both', 'emarker2', {idx_channel, 'o', 'r', 3, 2});
title(['oth hi-lo'])

subplot(3, 3, 7)
topoplot(topovals(7, :) - topovals(1, :), chanlocs, 'maplimits', plotlims, 'electrodes', 'on', 'style', 'both', 'emarker2', {idx_channel, 'o', 'r', 3, 2});
title(['oth-slf in lo'])

subplot(3, 3, 8)
topoplot(topovals(8, :) - topovals(2, :), chanlocs, 'maplimits', plotlims, 'electrodes', 'on', 'style', 'both', 'emarker2', {idx_channel, 'o', 'r', 3, 2});
title(['oth-slf in hi'])


% Plot tf for conditions
figure()

subplot(3, 3, 1)
pd = squeeze(mean(squeeze(ersps_patch(:, 1, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'linecolor','none')
colormap('jet')
set(gca,'clim', [-2, 2], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
colorbar;
title(['slf lo'])

subplot(3, 3, 2)
pd = squeeze(mean(squeeze(ersps_patch(:, 2, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'linecolor','none')
colormap('jet')
set(gca,'clim', [-2, 2], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
colorbar;
title(['slf hi'])

subplot(3, 3, 4)
pd = squeeze(mean(squeeze(ersps_patch(:, 7, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'linecolor','none')
colormap('jet')
set(gca,'clim', [-2, 2], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
colorbar;
title(['oth lo'])

subplot(3, 3, 5)
pd = squeeze(mean(squeeze(ersps_patch(:, 8, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'linecolor','none')
colormap('jet')
set(gca,'clim', [-2, 2], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
colorbar;
title(['oth hi'])

subplot(3, 3, 3)
pd = squeeze(mean(squeeze(ersps_patch(:, 2, :, :)), 1)) - squeeze(mean(squeeze(ersps_patch(:, 1, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'linecolor','none')
colormap('jet')
set(gca,'clim', [-2, 2], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
colorbar;
title(['slf hi-lo'])

subplot(3, 3, 6)
pd = squeeze(mean(squeeze(ersps_patch(:, 8, :, :)), 1)) - squeeze(mean(squeeze(ersps_patch(:, 7, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'linecolor','none')
colormap('jet')
set(gca,'clim', [-2, 2], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
colorbar;
title(['oth hi-lo'])

subplot(3, 3, 7)
pd = squeeze(mean(squeeze(ersps_patch(:, 7, :, :)), 1)) - squeeze(mean(squeeze(ersps_patch(:, 1, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'linecolor','none')
colormap('jet')
set(gca,'clim', [-2, 2], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
colorbar;
title(['oth-slf in lo'])

subplot(3, 3, 8)
pd = squeeze(mean(squeeze(ersps_patch(:, 8, :, :)), 1)) - squeeze(mean(squeeze(ersps_patch(:, 2, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 50, 'linecolor','none')
colormap('jet')
set(gca,'clim', [-2, 2], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
colorbar;
title(['oth-slf in hi'])

% Create output matrix and save (long format)
out_long = [];
counter = 0;

for s = 1 : size(ersps_vals, 1)
    for cond = 1 : 4
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
            val = ersps_vals(s, 7);
            fb = 2; % oth
            rew = 1; % lo
        elseif cond == 4
            val = ersps_vals(s, 8);
            fb = 2; % oth
            rew = 2; % hi
        end
        out_long(counter, :) = [s, fb, rew, val];
    end
end
fn = [PATH_TF_RESULTS, analysis_label, '.csv'];
writematrix(out_long, fn);