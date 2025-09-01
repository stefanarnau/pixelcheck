clear all

% Paths 
PATH_EEGLAB = '/home/plkn/eeglab2025.0.0/';
PATH_IN = '/mnt/data_dump/pixelcheck/2_cleaned_resplocked/';
PATH_TF_DATA = '/mnt/data_dump/pixelcheck/3_tf_data_resplocked_1stout/';

% Init EEGlab
addpath(PATH_EEGLAB);
eeglab;

% Get list of files
file_list = dir(fullfile(PATH_IN, '*.set'));

% Lade info
EEG = pop_loadset('filename', file_list(1).name, 'filepath', PATH_IN, 'loadmode', 'info');

% Set complex Morlet wavelet parameters
n_frq = 30;
frqrange = [2, 30];
tfres_range = [600, 240];

% Set wavelet time
wtime = -2 : 1 / EEG.srate : 2;

% Determine fft frqs
hz = linspace(0, EEG.srate, length(wtime));

% Create wavelet frequencies and tapering Gaussian widths in temporal domain
tf_freqs = logspace(log10(frqrange(1)), log10(frqrange(2)), n_frq);
fwhmTs = logspace(log10(tfres_range(1)), log10(tfres_range(2)), n_frq);

% Init matrices for wavelets
cmw = zeros(length(tf_freqs), length(wtime));
cmwX = zeros(length(tf_freqs), length(wtime));
tlim = zeros(1, length(tf_freqs));

% These will contain the wavelet widths as full width at 
% half maximum in the temporal and spectral domain
obs_fwhmT = zeros(1, length(tf_freqs));
obs_fwhmF = zeros(1, length(tf_freqs));

% Create the wavelets
for frq = 1 : length(tf_freqs)

    % Create wavelet with tapering gaussian corresponding to desired width in temporal domain
    cmw(frq, :) = exp(2 * 1i * pi * tf_freqs(frq) .* wtime) .* exp((-4 * log(2) * wtime.^2) ./ (fwhmTs(frq) / 1000)^2);

    % Normalize wavelet
    cmw(frq, :) = cmw(frq, :) ./ max(cmw(frq, :));

    % Create normalized freq domain wavelet
    cmwX(frq, :) = fft(cmw(frq, :)) ./ max(fft(cmw(frq, :)));

    % Determine observed fwhmT
    midt = dsearchn(wtime', 0);
    cmw_amp = abs(cmw(frq, :)) ./ max(abs(cmw(frq, :))); % Normalize cmw amplitude
    obs_fwhmT(frq) = wtime(midt - 1 + dsearchn(cmw_amp(midt : end)', 0.5)) - wtime(dsearchn(cmw_amp(1 : midt)', 0.5));

    % Determine observed fwhmF
    idx = dsearchn(hz', tf_freqs(frq));
    cmwx_amp = abs(cmwX(frq, :)); 
    obs_fwhmF(frq) = hz(idx - 1 + dsearchn(cmwx_amp(idx : end)', 0.5) - dsearchn(cmwx_amp(1 : idx)', 0.5));

end

% Define time window of analysis
prune_times = [-500, 1000]; 
tf_times = EEG.times(dsearchn(EEG.times', prune_times(1)) : dsearchn(EEG.times', prune_times(2)));

% Loop subjects
for s = 1 : length(file_list)

    % Result struct
    chanlocs = EEG.chanlocs;
    ersps = zeros(12, EEG.nbchan, length(tf_freqs), length(tf_times));
    itpcs = zeros(12, EEG.nbchan, length(tf_freqs), length(tf_times));

    % Load data
    EEG = pop_loadset('filename', file_list(s).name, 'filepath', PATH_IN, 'loadmode', 'all');

    % Code trialnumber in block
    EEG.trialinfo.trial_in_block = mod(EEG.trialinfo.trial_nr, 120);
    EEG.trialinfo.trial_in_block(EEG.trialinfo.trial_in_block == 0) = 120;

    % Find first blocks of each condition
    first_1 = 0;
    first_2 = 0;
    first_3 = 0;
    for t = 1 : height(EEG.trialinfo)
        if EEG.trialinfo.ma_condition(t) == 1 & first_1 == 0
            first_1 = EEG.trialinfo.block_nr(t);
        end
        if EEG.trialinfo.ma_condition(t) == 2 & first_2 == 0
            first_2 = EEG.trialinfo.block_nr(t);
        end
        if EEG.trialinfo.ma_condition(t) == 3 & first_3 == 0
            first_3 = EEG.trialinfo.block_nr(t);
        end
    end

    % Set trial exclusion criteria
    idx_keep = ~ismember(EEG.trialinfo.block_nr, [first_1, first_2, first_3]);

    % Exclude trials
    eeg_data = EEG.data(:, :, idx_keep);
    trialinfo = EEG.trialinfo(idx_keep, :);
    EEG.trials = sum(idx_keep);

    % Condition labels
    cond_label = {'slf err lo', 'slf err hi', 'slf flip lo', 'slf flip hi', 'slf corr lo', 'slf corr hi',...
                    'oth err lo', 'oth err hi', 'oth flip lo', 'oth flip hi', 'oth corr lo', 'oth corr hi'};

    % Loop conditions
    cond_idx = {};
    for cond = 1 : 12
        switch cond
            case 1 
                cond_idx{cond} = trialinfo.ma_condition == 2 & trialinfo.fb_correct == 0 & trialinfo.fb_flipped == 0 & trialinfo.reward_condition == 0;
            case 2 
                cond_idx{cond} = trialinfo.ma_condition == 2 & trialinfo.fb_correct == 0 & trialinfo.fb_flipped == 0 & trialinfo.reward_condition == 1;
            case 3
                cond_idx{cond} = trialinfo.ma_condition == 2 & trialinfo.fb_correct == 0 & trialinfo.fb_flipped == 1 & trialinfo.reward_condition == 0;
            case 4
                cond_idx{cond} = trialinfo.ma_condition == 2 & trialinfo.fb_correct == 0 & trialinfo.fb_flipped == 1 & trialinfo.reward_condition == 1;
            case 5
                cond_idx{cond} = trialinfo.ma_condition == 2 & trialinfo.fb_correct == 1 & trialinfo.reward_condition == 0;
            case 6
                cond_idx{cond} = trialinfo.ma_condition == 2 & trialinfo.fb_correct == 1 & trialinfo.reward_condition == 1;
            case 7
                cond_idx{cond} = trialinfo.ma_condition == 3 & trialinfo.fb_correct == 0 & trialinfo.fb_flipped == 0 & trialinfo.reward_condition == 0;
            case 8
                cond_idx{cond} = trialinfo.ma_condition == 3 & trialinfo.fb_correct == 0 & trialinfo.fb_flipped == 0 & trialinfo.reward_condition == 1;
            case 9
                cond_idx{cond} = trialinfo.ma_condition == 3 & trialinfo.fb_correct == 0 & trialinfo.fb_flipped == 1 & trialinfo.reward_condition == 0;
            case 10
                cond_idx{cond} = trialinfo.ma_condition == 3 & trialinfo.fb_correct == 0 & trialinfo.fb_flipped == 1 & trialinfo.reward_condition == 1;
            case 11
                cond_idx{cond} = trialinfo.ma_condition == 3 & trialinfo.fb_correct == 1 & trialinfo.reward_condition == 0;
            case 12
                cond_idx{cond} = trialinfo.ma_condition == 3 & trialinfo.fb_correct == 1 & trialinfo.reward_condition == 1;
        end

    end

    % Loop channels
    for ch = 1 : EEG.nbchan

        % Init tf matrices
        powcube = NaN(length(tf_freqs), EEG.pnts, EEG.trials);
        phacube = NaN(length(tf_freqs), EEG.pnts, EEG.trials);

        % Talk
        fprintf('\ntf-decomposition | subject %i/%i | channel %i/%i\n', s, length(file_list), ch, EEG.nbchan);

        % Get channel data
        channel_data = squeeze(eeg_data(ch, :, :));

        % convolution length
        convlen = size(channel_data, 1) * size(channel_data, 2) + size(cmw, 2) - 1;

        % cmw to freq domain and scale
        cmwX = zeros(length(tf_freqs), convlen);
        for f = 1 : length(tf_freqs)
            cmwX(f, :) = fft(cmw(f, :), convlen);
            cmwX(f, :) = cmwX(f, :) ./ max(cmwX(f, :));
        end

        % Get TF-power
        tmp = fft(reshape(channel_data, 1, []), convlen);
        for f = 1 : length(tf_freqs)
            as = ifft(cmwX(f, :) .* tmp); 
            as = as(((size(cmw, 2) - 1) / 2) + 1 : end - ((size(cmw, 2) - 1) / 2));
            as = reshape(as, EEG.pnts, EEG.trials);
            powcube(f, :, :) = abs(as) .^ 2;
            phacube(f, :, :) = exp(1i*angle(as)); 
        end
        
        % Cut edges
        powcube = powcube(:, dsearchn(EEG.times', prune_times(1)) : dsearchn(EEG.times', prune_times(2)), :);
        phacube = phacube(:, dsearchn(EEG.times', prune_times(1)) : dsearchn(EEG.times', prune_times(2)), :);

        % Get condition general baseline values
        ersp_bl = [-500, -200];
        tmp = squeeze(mean(powcube, 3));
        [~, blidx1] = min(abs(tf_times - ersp_bl(1)));
        [~, blidx2] = min(abs(tf_times - ersp_bl(2)));
        blvals = squeeze(mean(tmp(:, blidx1 : blidx2), 2));

        % Loop conditions and calculate ersp and itpc
        for cond = 1 : 12
            ersps(cond, ch, :, :) = single(10 * log10(bsxfun(@rdivide, squeeze(mean(powcube(:, :, cond_idx{cond}), 3)), blvals)));
            itpcs(cond, ch, :, :) = single(abs(mean(phacube(:, :, cond_idx{cond}), 3)));
        end

    end % end channel loop

    % Save
    save([PATH_TF_DATA, 'chanlocs.mat'], 'chanlocs');
    save([PATH_TF_DATA, 'tf_freqs.mat'], 'tf_freqs');
    save([PATH_TF_DATA, 'tf_times.mat'], 'tf_times');
    save([PATH_TF_DATA, num2str(trialinfo.id(1)), '_ersps.mat'], 'ersps');
    save([PATH_TF_DATA, num2str(trialinfo.id(1)), '_itpcs.mat'], 'itpcs');

end % end subject loop
