clear all;

% File path
PATH_IN = '/mnt/data_dump/pixelcheck_resplock_ged/feedback_locked/';
PATH_EEGLAB = '/home/plkn/eeglab2025.0.0/';
PATH_OUT = '/mnt/data_dump/pixelcheck_resplock_ged/ged_calculated_flipped_theta/';

% Init EEGlab
addpath(PATH_EEGLAB);
eeglab;

% Get list of files
file_list = dir(fullfile(PATH_IN, '*.set'));

res = [];

% Iterate files
for f = 1 : length(file_list)

    % Get filename
    fn = file_list(f).name;

    % Subject id
    subject = fn(1 : 5);

    % Load data
    EEG = pop_loadset('filename', fn, 'filepath', PATH_IN, 'loadmode', 'info');

    res(f, :) = [EEG.trials, size(EEG.trialinfo, 1), EEG.trials - size(EEG.trialinfo, 1)];

end





