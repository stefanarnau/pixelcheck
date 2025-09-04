clear all

% Paths 
PATH_EEGLAB = '/home/plkn/eeglab2025.0.0/';
PATH_FIELDTRIP = '/home/plkn/fieldtrip-master/';
PATH_TF_DATA = '/mnt/data_dump/pixelcheck/3_tf_data_cuelocked/';
PATH_TF_RESULTS = '/mnt/data_dump/pixelcheck/4_tf_results_cuelocked/';

% Init EEGlab
addpath(PATH_EEGLAB);
eeglab;

% Init fieldtrip
addpath(PATH_FIELDTRIP);
ft_defaults;

% Get list of files
file_list = dir(fullfile(PATH_TF_DATA, '*ersps.mat'));
  
% Get number of subjects
n_subjects = length(file_list);

% Loop subjects and load data
for s = 1 : n_subjects

    % Load metadata
    load([PATH_TF_DATA, 'chanlocs.mat']);
    load([PATH_TF_DATA, 'tf_freqs.mat']);
    load([PATH_TF_DATA, 'tf_times.mat']);

    % Build elec struct
    for ch = 1 : length(chanlocs)
        elec.label{ch} = chanlocs(ch).labels;
        elec.elecpos(ch, :) = [chanlocs(ch).X, chanlocs(ch).Y, chanlocs(ch).Z];
        elec.chanpos(ch, :) = [chanlocs(ch).X, chanlocs(ch).Y, chanlocs(ch).Z];
    end

    % Get id as string
    id_string = regexp(file_list(s).name, '\d+', 'match');
    id_string = id_string{1};

    % Load data
    load([PATH_TF_DATA, id_string, '_ersps.mat']); % cond x chan x freq x time

    % Prune to cti
    time_idx = tf_times >= -500 & tf_times <= 1800;
    tf_times = tf_times(time_idx);
    ersps = ersps(:, :, :, time_idx);

    % Main effect agency
    ersp_neu = squeeze(mean(ersps([1, 2], :, :, :), 1));
    ersp_slf = squeeze(mean(ersps([3, 4], :, :, :), 1));
    ersp_oth = squeeze(mean(ersps([5, 6], :, :, :), 1));

    % Main effect reward
    ersp_lo = squeeze(mean(ersps([1, 3, 5], :, :, :), 1));
    ersp_hi = squeeze(mean(ersps([2, 4, 6], :, :, :), 1));

    % Interaction
    ersp_int_neu = squeeze(ersps(2, :, :, :)) - squeeze(ersps(1, :, :, :));
    ersp_int_slf = squeeze(ersps(4, :, :, :)) - squeeze(ersps(3, :, :, :));
    ersp_int_oth = squeeze(ersps(6, :, :, :)) - squeeze(ersps(5, :, :, :));

    % Build fieldtrip structs 
    tf_ersp_neu.powspctrm = ersp_neu;
    tf_ersp_neu.dimord    = 'chan_freq_time';
    tf_ersp_neu.label     = elec.label;
    tf_ersp_neu.freq      = tf_freqs;
    tf_ersp_neu.time      = tf_times;

    tf_ersp_slf.powspctrm = ersp_slf;
    tf_ersp_slf.dimord    = 'chan_freq_time';
    tf_ersp_slf.label     = elec.label;
    tf_ersp_slf.freq      = tf_freqs;
    tf_ersp_slf.time      = tf_times;

    tf_ersp_oth.powspctrm = ersp_oth;
    tf_ersp_oth.dimord    = 'chan_freq_time';
    tf_ersp_oth.label     = elec.label;
    tf_ersp_oth.freq      = tf_freqs;
    tf_ersp_oth.time      = tf_times;

    tf_ersp_lo.powspctrm = ersp_lo;
    tf_ersp_lo.dimord    = 'chan_freq_time';
    tf_ersp_lo.label     = elec.label;
    tf_ersp_lo.freq      = tf_freqs;
    tf_ersp_lo.time      = tf_times;

    tf_ersp_hi.powspctrm = ersp_hi;
    tf_ersp_hi.dimord    = 'chan_freq_time';
    tf_ersp_hi.label     = elec.label;
    tf_ersp_hi.freq      = tf_freqs;
    tf_ersp_hi.time      = tf_times;

    tf_ersp_int_neu.powspctrm = ersp_int_neu;
    tf_ersp_int_neu.dimord    = 'chan_freq_time';
    tf_ersp_int_neu.label     = elec.label;
    tf_ersp_int_neu.freq      = tf_freqs;
    tf_ersp_int_neu.time      = tf_times;

    tf_ersp_int_slf.powspctrm = ersp_int_slf;
    tf_ersp_int_slf.dimord    = 'chan_freq_time';
    tf_ersp_int_slf.label     = elec.label;
    tf_ersp_int_slf.freq      = tf_freqs;
    tf_ersp_int_slf.time      = tf_times;

    tf_ersp_int_oth.powspctrm = ersp_int_oth;
    tf_ersp_int_oth.dimord    = 'chan_freq_time';
    tf_ersp_int_oth.label     = elec.label;
    tf_ersp_int_oth.freq      = tf_freqs;
    tf_ersp_int_oth.time      = tf_times;

    % Collect
    all_tf_ersp_neu{s} = tf_ersp_neu;
    all_tf_ersp_slf{s} = tf_ersp_slf;
    all_tf_ersp_oth{s} = tf_ersp_oth;
    all_tf_ersp_lo{s} = tf_ersp_lo;
    all_tf_ersp_hi{s} = tf_ersp_hi;
    all_tf_ersp_int_neu{s} = tf_ersp_int_neu;
    all_tf_ersp_int_slf{s} = tf_ersp_int_slf;
    all_tf_ersp_int_oth{s} = tf_ersp_int_oth;

end

% Calculate grand averages
cfg = [];
cfg.keepindividual = 'yes';
GA_ersp_neu = ft_freqgrandaverage(cfg, all_tf_ersp_neu{1, :});
GA_ersp_slf = ft_freqgrandaverage(cfg, all_tf_ersp_slf{1, :});
GA_ersp_oth = ft_freqgrandaverage(cfg, all_tf_ersp_oth{1, :});

GA_ersp_lo = ft_freqgrandaverage(cfg, all_tf_ersp_lo{1, :});
GA_ersp_hi = ft_freqgrandaverage(cfg, all_tf_ersp_hi{1, :});

GA_ersp_int_neu = ft_freqgrandaverage(cfg, all_tf_ersp_int_neu{1, :});
GA_ersp_int_slf = ft_freqgrandaverage(cfg, all_tf_ersp_int_slf{1, :});
GA_ersp_int_oth = ft_freqgrandaverage(cfg, all_tf_ersp_int_oth{1, :});

aa = bb;

% Prepare layout
cfg      = [];
cfg.elec = elec;
cfg.rotate = 90;
layout = ft_prepare_layout(cfg);

% Define neighbours
cfg                 = [];
cfg.layout          = layout;
cfg.feedback        = 'no';
cfg.method          = 'triangulation'; 
cfg.neighbours      = ft_prepare_neighbours(cfg, GA_ersp_neu);
neighbours          = cfg.neighbours;

% Testparams
testalpha   = 0.05;
voxelalpha  = 0.05;
nperm       = 1000;

% Omnibs ANOVAs
execute = 0;
if execute == 1

    % Set config for 2 level factor
    cfg = [];
    cfg.tail             = 0;
    cfg.statistic        = 'depsamplesT';
    cfg.alpha            = testalpha;
    cfg.neighbours       = neighbours;
    cfg.minnbchan        = 2;
    cfg.method           = 'montecarlo';
    cfg.correctm         = 'cluster';
    cfg.clustertail      = 0;
    cfg.clusteralpha     = voxelalpha;
    cfg.clusterstatistic = 'maxsum';
    cfg.numrandomization = nperm;
    cfg.computecritval   = 'yes'; 
    cfg.ivar             = 1;
    cfg.uvar             = 2;
    cfg.design   = [ones(1, n_subjects), 2 * ones(1, n_subjects); 1 : n_subjects, 1 : n_subjects];

    % Test main effect reward
    [stat_reward] = ft_freqstatistics(cfg, GA_ersp_lo, GA_ersp_hi);
    save([PATH_TF_RESULTS, 'stat_reward.mat'], 'stat_reward');

    % Set config for 3 level factor
    cfg = [];
    cfg.tail             = 1;
    cfg.statistic        = 'depsamplesFmultivariate';
    cfg.alpha            = testalpha;
    cfg.neighbours       = neighbours;
    cfg.minnbchan        = 2;
    cfg.method           = 'montecarlo';
    cfg.correctm         = 'cluster';
    cfg.clustertail      = 1;
    cfg.clusteralpha     = voxelalpha;
    cfg.clusterstatistic = 'maxsum';
    cfg.numrandomization = nperm;
    cfg.computecritval   = 'yes'; 
    cfg.ivar             = 1;
    cfg.uvar             = 2;
    cfg.design   = [ones(1, n_subjects), 2 * ones(1, n_subjects), 3 * ones(1, n_subjects); 1 : n_subjects, 1 : n_subjects, 1 : n_subjects];

    % Test main effect agency 
    [stat_agency] = ft_freqstatistics(cfg, GA_ersp_neu, GA_ersp_slf, GA_ersp_oth);
    save([PATH_TF_RESULTS, 'stat_agency.mat'], 'stat_agency');

    % Test interaction 
    [stat_interaction] = ft_freqstatistics(cfg, GA_ersp_int_neu, GA_ersp_int_slf, GA_ersp_int_oth);
    save([PATH_TF_RESULTS, 'stat_interaction.mat'], 'stat_interaction');
else

    % Re-load stat structs
    load([PATH_TF_RESULTS, 'stat_reward.mat']);
    load([PATH_TF_RESULTS, 'stat_agency.mat']);
    load([PATH_TF_RESULTS, 'stat_interaction.mat']);

end


% Set config for 2 level factor
cfg = [];
cfg.tail             = 0;
cfg.statistic        = 'depsamplesT';
cfg.alpha            = testalpha;
cfg.neighbours       = neighbours;
cfg.minnbchan        = 2;
cfg.method           = 'montecarlo';
cfg.correctm         = 'cluster';
cfg.clustertail      = 0;
cfg.clusteralpha     = voxelalpha;
cfg.clusterstatistic = 'maxsum';
cfg.numrandomization = nperm;
cfg.computecritval   = 'yes'; 
cfg.ivar             = 1;
cfg.uvar             = 2;
cfg.design   = [ones(1, n_subjects), 2 * ones(1, n_subjects); 1 : n_subjects, 1 : n_subjects];

% Test main effect reward
[stat_neu_v_slf] = ft_freqstatistics(cfg, GA_ersp_neu, GA_ersp_slf);
[stat_neu_v_oth] = ft_freqstatistics(cfg, GA_ersp_neu, GA_ersp_oth);
[stat_slf_v_oth] = ft_freqstatistics(cfg, GA_ersp_slf, GA_ersp_oth);
save([PATH_TF_RESULTS, 'stat_neu_v_slf.mat'], 'stat_neu_v_slf');
save([PATH_TF_RESULTS, 'stat_neu_v_oth.mat'], 'stat_neu_v_oth');
save([PATH_TF_RESULTS, 'stat_slf_v_oth.mat'], 'stat_slf_v_oth');









idx_cluster1 = stat_agency.posclusterslabelmat == 1;

idx_channel = sum(idx_cluster1, [2, 3]);
tf_mask_cluster1 = logical(squeeze(sum(idx_cluster1, 1)));

pow_neu = squeeze(mean(GA_ersp_neu.powspctrm, [1, 2]));
pow_slf = squeeze(mean(GA_ersp_slf.powspctrm, [1, 2]));
pow_oth = squeeze(mean(GA_ersp_oth.powspctrm, [1, 2]));



figure()
subplot(2, 2, 1)
contourf(tf_times, tf_freqs, pow_neu, 40, 'linecolor','none')
hold on
contour(tf_times, tf_freqs, tf_mask_cluster1, 1, 'linecolor', 'k', 'LineWidth', 2)
colormap('jet')
set(gca,'clim', [-2, 2], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
colorbar;
title(['neu'])

subplot(2, 2, 2)
contourf(tf_times, tf_freqs, pow_slf, 40, 'linecolor','none')
hold on
contour(tf_times, tf_freqs, tf_mask_cluster1, 1, 'linecolor', 'k', 'LineWidth', 2)
colormap('jet')
set(gca,'clim', [-2, 2], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
colorbar;
title(['slf'])

subplot(2, 2, 3)
contourf(tf_times, tf_freqs, pow_oth, 40, 'linecolor','none')
hold on
contour(tf_times, tf_freqs, tf_mask_cluster1, 1, 'linecolor', 'k', 'LineWidth', 2)
colormap('jet')
set(gca,'clim', [-2, 2], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
colorbar;
title(['oth'])

subplot(2, 2, 4)

pd_neu = mean(pow_neu(tf_freqs >= 8 & tf_freqs <= 12, :), 1);
pd_slf = mean(pow_slf(tf_freqs >= 8 & tf_freqs <= 12, :), 1);
pd_oth = mean(pow_oth(tf_freqs >= 8 & tf_freqs <= 12, :), 1);

plot(tf_times, [pd_neu; pd_slf; pd_oth]);



clust_thresh = 0.05;  % Schwellenwert für signifikante Cluster
clusts = struct();
cnt = 0;

stat = stat_ersp_macond;  % Hier ggf. stat_itpc_macond verwenden, falls nötig

if isfield(stat, 'posclusters') && ~isempty(stat.posclusters)
    pos_idx = find([stat.posclusters(1, :).prob] < clust_thresh);
    for c = 1 : numel(pos_idx)
        cnt = cnt + 1;
        clusts(cnt).testlabel = 'stat_ersp_macond';
        clusts(cnt).clustnum = pos_idx(c);
        clusts(cnt).time = stat.time;
        clusts(cnt).freq = stat.freq;
        clusts(cnt).prob = stat.posclusters(1, pos_idx(c)).prob;
        clusts(cnt).idx = stat.posclusterslabelmat == pos_idx(c);  % Maske für Cluster
        clusts(cnt).stats = stat.stat;
        clusts(cnt).chans_sig = find(logical(mean(clusts(cnt).idx, [2, 3])));
    end
else
    warning('Keine signifikanten Cluster in stat_ersp_macond gefunden.');
end

    % Calculate adjusted partial eta squared
    apes_macond = [];
    for ch = 1 : 65
        petasq = (squeeze(stat_ersp_macond.stat(ch, :, :)) .* 2) ./ ((squeeze(stat_ersp_macond .stat(ch, :, :)) * 2) + (n_subjects - 1));
        adj_petasq = petasq - (1 - petasq) .* (2 / (n_subjects - 1));
        apes_macond(ch, :, :) = adj_petasq;
    end

    % Plot ersps for ma condition
    idx_channel = logical(squeeze(mean(stat_ersp_macond.posclusterslabelmat == 1, [2, 3])));

    figure()
    subplot(2, 2, 1)
    pd = squeeze(mean(squeeze(mean(squeeze(itpcs(1, :, idx_channel, :, :)), 1)), 1));
    contourf(tf_times, tf_freqs, pd, 40, 'linecolor','none')
    hold on
    contour(tf_times, tf_freqs, logical(squeeze(mean(clusts(cnt).idx, 1))), 1, 'linecolor', 'k', 'LineWidth', 2)
    colormap('jet')
    set(gca,'clim', [-1, 1], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
    colorbar;
    title(['neu'])

    subplot(2, 2, 2)
    pd = squeeze(mean(squeeze(mean(squeeze(itpcs(2, :, idx_channel, :, :)), 1)), 1));
    contourf(tf_times, tf_freqs, pd, 40, 'linecolor','none')
    hold on
    contour(tf_times, tf_freqs, logical(squeeze(mean(clusts(cnt).idx, 1))), 1, 'linecolor', 'k', 'LineWidth', 2)
    colormap('jet')
    set(gca,'clim', [-1, 1], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
    colorbar;
    title(['slf'])

    subplot(2, 2, 3)
    pd = squeeze(mean(squeeze(mean(squeeze(itpcs(3, :, idx_channel, :, :)), 1)), 1));
    contourf(tf_times, tf_freqs, pd, 40, 'linecolor','none')
    hold on
    contour(tf_times, tf_freqs, logical(squeeze(mean(clusts(cnt).idx, 1))), 1, 'linecolor', 'k', 'LineWidth', 2)
    colormap('jet')
    set(gca,'clim', [-1, 1], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
    colorbar;
    title(['oth'])

    subplot(2, 2, 4)
    pd = squeeze(mean(apes_macond(idx_channel, :, :), 1));
    contourf(tf_times, tf_freqs, pd, 40, 'linecolor','none')
    hold on
    contour(tf_times, tf_freqs, logical(squeeze(mean(clusts(cnt).idx, 1))), 1, 'linecolor', 'k', 'LineWidth', 2)
    colormap('jet')
    set(gca,'clim', [-0.8, 0.8], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
    colorbar;
    title(['effect sizes'])


    % ===========================================================================================================

    % Identify significant clusters
    clust_thresh = 0.05;
    clusts = struct();
    cnt = 0;
    stat_names = {'stat_ersp_macond'};
    for s = 1 : numel(stat_names)
        stat = eval(stat_names{s});
        if ~isempty(stat.posclusters)
            pos_idx = find([stat.posclusters(1, :).prob] < clust_thresh);
            for c = 1 : numel(pos_idx)
                cnt = cnt + 1;
                clusts(cnt).testlabel = stat_names{s};
                clusts(cnt).clustnum = cnt;
                clusts(cnt).time = stat.time;
                clusts(cnt).freq = stat.freq;
                clusts(cnt).prob = stat.posclusters(1, pos_idx(c)).prob;
                clusts(cnt).idx = stat.posclusterslabelmat == pos_idx(c);
                clusts(cnt).stats = stat.stat;
                clusts(cnt).chans_sig = find(logical(mean(clusts(cnt).idx, [2, 3])));
            end
        end
    end

    % Plot identified cluster
    clinecol = 'k';
    cmap = 'jet';
    for cnt = 1 : numel(clusts)

        figure('Visible', 'off'); clf;

        subplot(2, 2, 1)
        pd = squeeze(sum(clusts(cnt).stats, 1));
        contourf(clusts(cnt).time, clusts(cnt).freq, pd, 40, 'linecolor','none')
        hold on
        contour(clusts(cnt).time, clusts(cnt).freq, logical(squeeze(mean(clusts(cnt).idx, 1))), 1, 'linecolor', clinecol, 'LineWidth', 2)
        colormap(cmap)
        set(gca, 'xlim', [clusts(cnt).time(1), clusts(cnt).time(end)], 'clim', [-max(abs(pd(:))), max(abs(pd(:)))], 'YScale', 'lin', 'YTick', [4, 8, 12, 20])
        colorbar;
        title(['sum t across chans'], 'FontSize', 10)

        subplot(2, 2, 2)
        pd = squeeze(mean(clusts(cnt).idx, 1));
        contourf(clusts(cnt).time, clusts(cnt).freq, pd, 40, 'linecolor','none')
        hold on
        contour(clusts(cnt).time, clusts(cnt).freq, logical(squeeze(mean(clusts(cnt).idx, 1))), 1, 'linecolor', clinecol, 'LineWidth', 2)
        colormap(cmap)
        set(gca, 'xlim', [clusts(cnt).time(1), clusts(cnt).time(end)], 'clim', [-1, 1], 'YScale', 'lin', 'YTick', [4, 8, 12, 20])
        colorbar;
        title(['proportion chans significant'], 'FontSize', 10)

        subplot(2, 2, 3)
        pd = squeeze(sum(clusts(cnt).stats, [2, 3]));
        topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'on');
        colormap(cmap)
        set(gca, 'clim', [-max(abs(pd(:))), max(abs(pd(:)))])
        colorbar;
        title(['sum t per electrode'], 'FontSize', 10)

        subplot(2, 2, 4)
        pd = squeeze(mean(clusts(cnt).idx, [2, 3]));
        topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'on');
        colormap(cmap)
        set(gca, 'clim', [-1, 1])
        colorbar;
        title(['proportion tf-points significant'], 'FontSize', 10)

        saveas(gcf, [PATH_FIRST_PLOTS 'clustnum_' num2str(clusts(cnt).clustnum) '_' clusts(cnt).testlabel '.png']); 
    end



