clear all

% Paths 
PATH_EEGLAB = '/home/plkn/eeglab2025.0.0/';
PATH_AUTOCLEANED_ERP = '/mnt/data_dump/pixelcheck/2_cleaned_resplocked_erp/';
PATH_RESPPREP = '/mnt/data_dump/pixelcheck/3_respprep/';

% Subject list 
subject_list = {'VP001', 'VP003', 'VP004', 'VP005', 'VP006', 'VP007', 'VP009', 'VP010', 'VP011', ...)
                'VP012', 'VP013', 'VP014', 'VP015', 'VP016', 'VP017', 'VP018', 'VP019', 'VP020', 'VP021', ...
                'VP022', 'VP023', 'VP024', 'VP025', 'VP026', 'VP027', 'VP028', 'VP029', 'VP030', 'VP0301', ...
                'VP032', 'VP033_1', 'VP034', 'VP035', 'VP036', 'VP037', 'VP038', 'VP041'};

% Initialize EEGLab
addpath(PATH_EEGLAB);
eeglab;

% An erp container
all_erp = [];
all_ersp_raw = [];
all_ersp_specific = [];
all_ersp_general = [];
all_itpc = [];

% Iterate subjects
n_trial = [];
for s = 1 : length(subject_list)

    % Get id
    id = str2double(subject_list{s}(3 : 5));

    % Load data
    EEG = pop_loadset('filename', [num2str(id), '_cleaned.set'], 'filepath', PATH_AUTOCLEANED_ERP, 'loadmode', 'all');

    EEG = pop_select(EEG, 'channel', [5, 6, 35, 36, 43, 44]);

    % correct_key_color:
    % 0 means even enum, means color2
    % 1 means odd enum, means color1

    % id:

    % odd means yellow left
    % even means blue left

    % color1 is always the left button, color 2 is always the right button
    % in trialinfo, color1 os coded as 1 and color2 is coded as 0

    % Loop conditions and collect indices of correct trials
    cond_idx = {};
    for cond = 1 : 12
        switch cond
            case 1 
                cond_idx{cond} = find(EEG.trialinfo.ma_condition == 1 & EEG.trialinfo.reward_condition == 0 & EEG.trialinfo.accuracy == 1 & EEG.trialinfo.correct_key_color == 1);
                n_trial(s, cond) = length(cond_idx{cond});
            case 2 
                cond_idx{cond} = find(EEG.trialinfo.ma_condition == 1 & EEG.trialinfo.reward_condition == 1 & EEG.trialinfo.accuracy == 1 & EEG.trialinfo.correct_key_color == 1);
                n_trial(s, cond) = length(cond_idx{cond});
            case 3
                cond_idx{cond} = find(EEG.trialinfo.ma_condition == 2 & EEG.trialinfo.reward_condition == 0 & EEG.trialinfo.accuracy == 1 & EEG.trialinfo.correct_key_color == 1);
                n_trial(s, cond) = length(cond_idx{cond});
            case 4
                cond_idx{cond} = find(EEG.trialinfo.ma_condition == 2 & EEG.trialinfo.reward_condition == 1 & EEG.trialinfo.accuracy == 1 & EEG.trialinfo.correct_key_color == 1);
                n_trial(s, cond) = length(cond_idx{cond});
            case 5
                cond_idx{cond} = find(EEG.trialinfo.ma_condition == 3 & EEG.trialinfo.reward_condition == 0 & EEG.trialinfo.accuracy == 1 & EEG.trialinfo.correct_key_color == 1);
                n_trial(s, cond) = length(cond_idx{cond});
            case 6
                cond_idx{cond} = find(EEG.trialinfo.ma_condition == 3 & EEG.trialinfo.reward_condition == 1 & EEG.trialinfo.accuracy == 1 & EEG.trialinfo.correct_key_color == 1);
                n_trial(s, cond) = length(cond_idx{cond});
            case 7 
                cond_idx{cond} = find(EEG.trialinfo.ma_condition == 1 & EEG.trialinfo.reward_condition == 0 & EEG.trialinfo.accuracy == 1 & EEG.trialinfo.correct_key_color == 0);
                n_trial(s, cond) = length(cond_idx{cond});
            case 8 
                cond_idx{cond} = find(EEG.trialinfo.ma_condition == 1 & EEG.trialinfo.reward_condition == 1 & EEG.trialinfo.accuracy == 1 & EEG.trialinfo.correct_key_color == 0);
                n_trial(s, cond) = length(cond_idx{cond});
            case 9
                cond_idx{cond} = find(EEG.trialinfo.ma_condition == 2 & EEG.trialinfo.reward_condition == 0 & EEG.trialinfo.accuracy == 1 & EEG.trialinfo.correct_key_color == 0);
                n_trial(s, cond) = length(cond_idx{cond});
            case 10
                cond_idx{cond} = find(EEG.trialinfo.ma_condition == 2 & EEG.trialinfo.reward_condition == 1 & EEG.trialinfo.accuracy == 1 & EEG.trialinfo.correct_key_color == 0);
                n_trial(s, cond) = length(cond_idx{cond});
            case 11
                cond_idx{cond} = find(EEG.trialinfo.ma_condition == 3 & EEG.trialinfo.reward_condition == 0 & EEG.trialinfo.accuracy == 1 & EEG.trialinfo.correct_key_color == 0);
                n_trial(s, cond) = length(cond_idx{cond});
            case 12
                cond_idx{cond} = find(EEG.trialinfo.ma_condition == 3 & EEG.trialinfo.reward_condition == 1 & EEG.trialinfo.accuracy == 1 & EEG.trialinfo.correct_key_color == 0);
                n_trial(s, cond) = length(cond_idx{cond});
        end
    end

    % TF decomp
    [ersp_raw, ersp_db_general, ersp_db_specific, itpc, tf_time, tf_frqs] = get_event_related_tf(EEG, cond_idx);

    % Collect
    all_ersp_raw(s, :, :, :, :) = ersp_raw;
    all_ersp_specific(s, :, :, :, :) = ersp_db_specific;
    all_ersp_general(s, :, :, :, :) = ersp_db_general;
    all_itpc(s, :, :, :, :) = itpc;

end

% Save measures
save([PATH_RESPPREP, 'all_erp.mat'], "all_erp")
save([PATH_RESPPREP, 'all_ersp_raw.mat'], "all_ersp_raw")
save([PATH_RESPPREP, 'all_ersp_specific.mat'], "all_ersp_specific")
save([PATH_RESPPREP, 'all_ersp_general.mat'], "all_ersp_general")
save([PATH_RESPPREP, 'all_itpc.mat'], "all_itpc")

lat_indices = [];
for s = 1 : length(subject_list)
    for cond = 1 : 6

        ipsi_left   = squeeze(all_ersp_raw(s, cond, [1, 3, 5], :, :));
        contra_left = squeeze(all_ersp_raw(s, cond, [2, 4, 6], :, :));
        ipsi_right   = squeeze(all_ersp_raw(s, cond + 6, [2, 4, 6], :, :));
        contra_right = squeeze(all_ersp_raw(s, cond + 6, [1, 3, 5], :, :));

        ipsi = (ipsi_left + ipsi_right) ./ 2;
        contra = (contra_left + contra_right) ./ 2;

        lat_indices(s, cond, :, :, :) = (ipsi - contra) ./ (ipsi + contra);

    end
end

% Remove subject in position 24
lat_indices(24, :, :, :, :) = [];

idx_frq = tf_frqs > 20;
lat_indices_mu = squeeze(mean(lat_indices(:, :, :, idx_frq, :), 4));

figure()
titles = {'nlo', 'nhi', 'slo', 'shi', 'olo', 'ohi'}
for cond = 1 : 6
    subplot(3, 2, cond)
    pd = squeeze(mean(squeeze(lat_indices(:, cond, 3, :, :)), 1));
    contourf(tf_time, tf_frqs, pd, 50, 'linecolor','none')
    colormap('jet')
    set(gca,'clim', [-0.1, 0.1], 'YScale', 'log', 'YTick', [4, 8, 12, 20])
    title(titles{cond})
end

figure()
pd = squeeze(mean(lat_indices_mu(:, :, [1, 2, 3], :), [1, 3]));
plot(tf_time, pd)
legend({'nlo', 'nhi', 'slo', 'shi', 'olo', 'ohi'});