clear all

% Paths 
PATH_EEGLAB = '/home/plkn/eeglab2025.0.0/';
PATH_FIELDTRIP = '/home/plkn/fieldtrip-master/';
PATH_TF_DATA = '/mnt/data_dump/pixelcheck/3_tf_data_resplocked/';
PATH_TF_RESULTS = '/mnt/data_dump/pixelcheck/4_tf_results_resplocked/';

% Init EEGlab
addpath(PATH_EEGLAB);
eeglab;

% Init fieldtrip
addpath(PATH_FIELDTRIP);
ft_defaults;

% Get list of files
file_list = dir(fullfile(PATH_IN, '*ersps.mat'));

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
                
% Get number of subjects
n_subjects = length(file_list);

% Loop subjects
for s = 1 : n_subjects

    % Get participant id
    id_string = regexp(file_list(s).name, '\d+', 'match'); 

    % Load
    load([PATH_TF_DATA, id_string, '_ersps.mat']); % cond x chan x freq x time

    % Compute difference ersps for error monitoring
    em_slf_lo = squeeze(ersps(1, :, :, :)) - squeeze(ersps(5, :, :, :));
    em_slf_hi = squeeze(ersps(2, :, :, :)) - squeeze(ersps(6, :, :, :));
    em_oth_lo = squeeze(ersps(7, :, :, :)) - squeeze(ersps(11, :, :, :));
    em_oth_hi = squeeze(ersps(8, :, :, :)) - squeeze(ersps(12, :, :, :));

    % Compute difference ersps for feedback flip
    ff_slf_lo = squeeze(ersps(3, :, :, :)) - squeeze(ersps(5, :, :, :));
    ff_slf_hi = squeeze(ersps(4, :, :, :)) - squeeze(ersps(6, :, :, :));
    ff_oth_lo = squeeze(ersps(9, :, :, :)) - squeeze(ersps(11, :, :, :));
    ff_oth_hi = squeeze(ersps(10, :, :, :)) - squeeze(ersps(12, :, :, :));

    % Build fieldtrip structs for error monitoring
    tf_em_slf.powspctrm = (em_slf_lo + em_slf_hi) ./ 2;
    tf_em_slf.dimord    = 'chan_freq_time';
    tf_em_slf.label     = elec.label;
    tf_em_slf.freq      = tf_freqs;
    tf_em_slf.time      = tf_times;

    tf_em_oth.powspctrm = (em_oth_lo + em_oth_hi) ./ 2;
    tf_em_oth.dimord    = 'chan_freq_time';
    tf_em_oth.label     = elec.label;
    tf_em_oth.freq      = tf_freqs;
    tf_em_oth.time      = tf_times;

    tf_em_lo.powspctrm = (em_slf_lo + em_oth_lo) ./ 2;
    tf_em_lo.dimord    = 'chan_freq_time';
    tf_em_lo.label     = elec.label;
    tf_em_lo.freq      = tf_freqs;
    tf_em_lo.time      = tf_times;

    tf_em_hi.powspctrm = (em_slf_hi + em_oth_hi) ./ 2;
    tf_em_hi.dimord    = 'chan_freq_time';
    tf_em_hi.label     = elec.label;
    tf_em_hi.freq      = tf_freqs;
    tf_em_hi.time      = tf_times;

    tf_em_diff_slf.powspctrm = em_slf_hi - em_slf_lo;
    tf_em_diff_slf.dimord    = 'chan_freq_time';
    tf_em_diff_slf.label     = elec.label;
    tf_em_diff_slf.freq      = tf_freqs;
    tf_em_diff_slf.time      = tf_times;

    tf_em_diff_oth.powspctrm = em_oth_hi - em_oth_lo;
    tf_em_diff_oth.dimord    = 'chan_freq_time';
    tf_em_diff_oth.label     = elec.label;
    tf_em_diff_oth.freq      = tf_freqs;
    tf_em_diff_oth.time      = tf_times;

    % Collect
    all_tf_em_slf{s} = tf_em_slf;
    all_tf_em_oth{s} = tf_em_oth;
    all_tf_em_lo{s} = tf_em_lo;
    all_tf_em_hi{s} = tf_em_hi;
    all_tf_em_diff_slf{s} = tf_em_diff_slf;
    all_tf_em_diff_oth{s} = tf_em_diff_oth;

    % Build fieldtrip structs for feedback flip
    tf_ff_slf.powspctrm = (ff_slf_lo + ff_slf_hi) ./ 2;
    tf_ff_slf.dimord    = 'chan_freq_time';
    tf_ff_slf.label     = elec.label;
    tf_ff_slf.freq      = tf_freqs;
    tf_ff_slf.time      = tf_times;

    tf_ff_oth.powspctrm = (ff_oth_lo + ff_oth_hi) ./ 2;
    tf_ff_oth.dimord    = 'chan_freq_time';
    tf_ff_oth.label     = elec.label;
    tf_ff_oth.freq      = tf_freqs;
    tf_ff_oth.time      = tf_times;

    tf_ff_lo.powspctrm = (ff_slf_lo + ff_oth_lo) ./ 2;
    tf_ff_lo.dimord    = 'chan_freq_time';
    tf_ff_lo.label     = elec.label;
    tf_ff_lo.freq      = tf_freqs;
    tf_ff_lo.time      = tf_times;

    tf_ff_hi.powspctrm = (ff_slf_hi + ff_oth_hi) ./ 2;
    tf_ff_hi.dimord    = 'chan_freq_time';
    tf_ff_hi.label     = elec.label;
    tf_ff_hi.freq      = tf_freqs;
    tf_ff_hi.time      = tf_times;

    tf_ff_diff_slf.powspctrm = ff_slf_hi - ff_slf_lo;
    tf_ff_diff_slf.dimord    = 'chan_freq_time';
    tf_ff_diff_slf.label     = elec.label;
    tf_ff_diff_slf.freq      = tf_freqs;
    tf_ff_diff_slf.time      = tf_times;

    tf_ff_diff_oth.powspctrm = ff_oth_hi - ff_oth_lo;
    tf_ff_diff_oth.dimord    = 'chan_freq_time';
    tf_ff_diff_oth.label     = elec.label;
    tf_ff_diff_oth.freq      = tf_freqs;
    tf_ff_diff_oth.time      = tf_times;

    % Collect
    all_tf_ff_slf{s} = tf_ff_slf;
    all_tf_ff_oth{s} = tf_ff_oth;
    all_tf_ff_lo{s} = tf_ff_lo;
    all_tf_ff_hi{s} = tf_ff_hi;
    all_tf_ff_diff_slf{s} = tf_ff_diff_slf;
    all_tf_ff_diff_oth{s} = tf_ff_diff_oth;

end % end subject loop

% Calculate grand averages
cfg = [];
cfg.keepindividual = 'yes';
GA_tf_em_slf = ft_freqgrandaverage(cfg, all_tf_em_slf{1, :});
GA_tf_em_oth = ft_freqgrandaverage(cfg, all_tf_em_oth{1, :});
GA_tf_em_lo = ft_freqgrandaverage(cfg, all_tf_em_lo{1, :});
GA_tf_em_hi = ft_freqgrandaverage(cfg, all_tf_em_hi{1, :});
GA_tf_em_diff_slf = ft_freqgrandaverage(cfg, all_tf_em_diff_slf{1, :});
GA_tf_em_diff_oth = ft_freqgrandaverage(cfg, all_tf_em_diff_oth{1, :});
GA_tf_ff_slf = ft_freqgrandaverage(cfg, all_tf_ff_slf{1, :});
GA_tf_ff_oth = ft_freqgrandaverage(cfg, all_tf_ff_oth{1, :});
GA_tf_ff_lo = ft_freqgrandaverage(cfg, all_tf_ff_lo{1, :});
GA_tf_ff_hi = ft_freqgrandaverage(cfg, all_tf_ff_hi{1, :});
GA_tf_ff_diff_slf = ft_freqgrandaverage(cfg, all_tf_ff_diff_slf{1, :});
GA_tf_ff_diff_oth = ft_freqgrandaverage(cfg, all_tf_ff_diff_oth{1, :});

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
voxelalpha  = 0.01;
nperm       = 1000;

% Set config. Same for all tests
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
cfg.design   = [ones(1, n_subjects), 2 * ones(1, n_subjects); 1 : n_subjects, 1 : n_subjects];

% Test main effect agency for error monitoring
[stat_em_agency] = ft_freqstatistics(cfg, GA_tf_em_slf, GA_tf_em_oth);

% Test main effect reward for error monitoring
[stat_em_reward] = ft_freqstatistics(cfg, GA_tf_em_lo, GA_tf_em_hi);

% Test interaction for error monitoring
[stat_em_interaction] = ft_freqstatistics(cfg, GA_tf_em_diff_slf, GA_tf_em_diff_oth);

% Test main effect agency for feedback flip
[stat_ff_agency] = ft_freqstatistics(cfg, GA_tf_ff_slf, GA_tf_ff_oth);

% Test main effect reward for feedback flip
[stat_ff_reward] = ft_freqstatistics(cfg, GA_tf_ff_lo, GA_tf_ff_hi);

% Test interaction for feedback flip
[stat_ff_interaction] = ft_freqstatistics(cfg, GA_tf_ff_diff_slf, GA_tf_ff_diff_oth);

% Save stat structs
save([PATH_TF_RESULTS, 'stat_em_agency.mat'], 'stat_em_agency');
save([PATH_TF_RESULTS, 'stat_em_reward.mat'], 'stat_em_reward');
save([PATH_TF_RESULTS, 'stat_em_interaction.mat'], 'stat_em_interaction');
save([PATH_TF_RESULTS, 'stat_ff_agency.mat'], 'stat_ff_agency');
save([PATH_TF_RESULTS, 'stat_ff_reward.mat'], 'stat_ff_reward');
save([PATH_TF_RESULTS, 'stat_ff_interaction.mat'], 'stat_ff_interaction');

aa=bb

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



