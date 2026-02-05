clear all

% Paths 
PATH_EEGLAB = '/home/plkn/eeglab2025.0.0/';
PATH_AUTOCLEANED_TF_CUE = '/mnt/data_dump/pixelcheck/2_cleaned_cuelocked_tf/';
PATH_AUTOCLEANED_TF_RESP = '/mnt/data_dump/pixelcheck/2_cleaned_resplocked_tf/';

% Get list of files
file_list_cue = dir(fullfile(PATH_AUTOCLEANED_TF_CUE, '*.set'));
file_list_resp = dir(fullfile(PATH_AUTOCLEANED_TF_RESP, '*.set'));

% Initialize EEGLab
addpath(PATH_EEGLAB);
eeglab;

stats_cue = [];
stats_resp = [];

% Iterate subjects
for s = 1 : length(file_list_cue)

    % Load data
    EEG = pop_loadset('filename', file_list_cue(s).name, 'filepath', PATH_AUTOCLEANED_TF_CUE, 'loadmode', 'info');

    n_chans = numel(EEG.chans_rejected_combined);
    n_epochs = numel(EEG.rejected_epochs);
    n_ics = numel(EEG.nobrainer);

    stats_cue(s, :) = [n_chans, n_epochs, n_ics];

end

mean(stats_cue, 1)
std(stats_cue, [], 1)

% Iterate subjects
for s = 1 : length(file_list_resp)

    % Load data
    EEG = pop_loadset('filename', file_list_resp(s).name, 'filepath', PATH_AUTOCLEANED_TF_RESP, 'loadmode', 'info');

    n_chans = numel(EEG.chans_rejected_combined);
    n_epochs = numel(EEG.rejected_epochs);
    n_ics = numel(EEG.nobrainer);

    stats_resp(s, :) = [n_chans, n_epochs, n_ics];

end

mean(stats_resp, 1)
std(stats_resp, [], 1)