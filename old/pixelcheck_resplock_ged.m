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

% Iterate files
for f = 1 : length(file_list)

    % Get filename
    fn = file_list(f).name;

    % Subject id
    subject = fn(1 : 5);

    % Load data
    EEG = pop_loadset('filename', fn, 'filepath', PATH_IN, 'loadmode', 'all');

    % To double precision
    eeg_data = double(EEG.data);
                
    % Exclude trials
    % TODO

    % Construct filter
    nyquist = EEG.srate / 2;
    h_pass = 4;
    l_pass = 7;
    transition_width = 0.2;
    filter_order = round(3 * (EEG.srate / h_pass));
    filter_freqs = [0, (1 - transition_width) * h_pass, h_pass, l_pass, (1 + transition_width) * l_pass, nyquist] / nyquist; 
    filter_response = [0, 0, 1, 1, 0, 0];
    filter_weights = firls(filter_order, filter_freqs, filter_response);

    % Reshape to 2d
    eeg_data_2d = reshape(eeg_data, [EEG.nbchan, EEG.pnts * EEG.trials]);

    % Apply filter
    eeg_data_filtered_2d = zeros(size(eeg_data_2d));
    for ch = 1 : EEG.nbchan
        eeg_data_filtered_2d(ch, :) = filtfilt(filter_weights, 1, eeg_data(ch, :));
    end

    % Reshape back to 3d
    eeg_data_filtered = reshape(eeg_data_filtered_2d, [EEG.nbchan, EEG.pnts, EEG.trials]);

    % Crop in time to remove edge artifacts
    crop_idx = EEG.times >= -500 & EEG.times <= 1500;
    times = EEG.times(crop_idx);
    eeg_data_filtered = eeg_data_filtered(:, crop_idx, :);
    eeg_data = eeg_data(:, crop_idx, :);

    % Get trial indices
    idx_correct = find(EEG.trialinfo.ma_condition > 1 & EEG.trialinfo.fb_correct == 1);
    idx_incorrect = find(EEG.trialinfo.ma_condition > 1 & EEG.trialinfo.fb_correct == 0 & EEG.trialinfo.fb_flipped == 0);
    idx_flipped = find(EEG.trialinfo.ma_condition > 1 & EEG.trialinfo.fb_correct == 0 & EEG.trialinfo.fb_flipped == 1);

    % Find indices of time points 
    tidx = (times >= 150 & times <= 750);

    % Init arrays for trial-specific covariance matrices
    covmats_correct = zeros(length(idx_correct), size(eeg_data, 1), size(eeg_data, 1));
    covmats_incorrect = zeros(length(idx_incorrect), size(eeg_data, 1), size(eeg_data, 1));
    covmats_flipped = zeros(length(idx_flipped), size(eeg_data, 1), size(eeg_data, 1));

    % Split eeg data
    eeg_data_correct = eeg_data_filtered(:, :, idx_correct);
    eeg_data_incorrect = eeg_data_filtered(:, :, idx_incorrect);
    eeg_data_flipped = eeg_data_filtered(:, :, idx_flipped);

    % Covariance matrix for each trial
    for trial_idx = 1 : length(idx_correct)
        tmp = squeeze(eeg_data_correct(:, tidx, trial_idx));
        covmats_correct(trial_idx, :, :) = cov(tmp');
    end
    for trial_idx = 1 : length(idx_incorrect)
        tmp = squeeze(eeg_data_incorrect(:, tidx, trial_idx));
        covmats_incorrect(trial_idx, :, :) = cov(tmp');
    end
    for trial_idx = 1 : length(idx_flipped)
        tmp = squeeze(eeg_data_flipped(:, tidx, trial_idx));
        covmats_flipped(trial_idx, :, :) = cov(tmp');
    end

    % Compute average covariance matrices
    ave_covmat_correct = squeeze(mean(covmats_correct, 1));
    ave_covmat_incorrect = squeeze(mean(covmats_incorrect, 1));
    ave_covmat_flipped = squeeze(mean(covmats_flipped, 1));

    % Apply shrinkage regularization to reference matrices
    g = 0.2;
    a = mean(eig(ave_covmat_correct));
    ave_covmat_correct = (1 - g) * ave_covmat_correct + g * a * eye(EEG.nbchan);

    % GED 
    [evecs, evals] = eig(ave_covmat_flipped, ave_covmat_correct);

    % Get minimum number of trials
    min_n = min([size(covmats_correct, 1), size(covmats_flipped, 1)]);

    % Create H0-distribution of eigenvalues
    nperm = 1000;
    eval_dist = zeros(nperm, 1);
    for p = 1 : nperm

        fprintf('Permutation %i/%i\n', p, nperm);

        % Draw random samples
        sample_correct_1 = datasample(covmats_correct, floor(min_n / 2), 1);
        sample_correct_2 = datasample(covmats_correct, floor(min_n / 2), 1);
        sample_flipped_1 = datasample(covmats_flipped, floor(min_n / 2), 1);
        sample_flipped_2 = datasample(covmats_flipped, floor(min_n / 2), 1);

        % Get S and R covmats from permuted data
        random_correct = squeeze(mean(cat(1, sample_correct_1, sample_flipped_1), 1));
        random_flipped = squeeze(mean(cat(1, sample_correct_2, sample_flipped_2), 1));

        % Apply shrinkage regularization to reference matrices
        g = 0.1;
        a = mean(eig(random_correct));
        random_correct = (1 - g) * random_correct + g * a * eye(EEG.nbchan);

        % GED 
        [~, random_evals] = eig(random_flipped, random_correct);

        % Get largest eigenvalue
        eval_dist(p) = max(diag(random_evals));

    end

    % Sort
    [evals, srtidx] = sort(diag(evals), 'descend');
    evecs = evecs(:, srtidx);
    eval_dist = sort(eval_dist, 'descend');

    % Get 95 percentile of eigenvalue distribution
    p95_threshold = prctile(eval_dist, 95);

    % Normalize eigenvectors
    evecs = bsxfun(@rdivide, evecs, sqrt(sum(evecs .^ 2, 1)));

    % Iterate components
    ged_maps = zeros(EEG.nbchan, EEG.nbchan);
    ged_time_series = zeros(EEG.nbchan, length(EEG.times), EEG.trials);
    ged_time_series_filtered = zeros(EEG.nbchan, length(EEG.times), EEG.trials);

    for cmp = 1 : EEG.nbchan

        % Compute map
        ged_maps(cmp, :) = evecs(:, cmp)' * ave_covmat_flipped;

        % Flip map if necessary
        [~, idx] = max(abs(ged_maps(cmp, :)));
        ged_maps(cmp, :) = ged_maps(cmp, :) * sign(ged_maps(cmp, idx));

        % Compute time series for component
        ged_time_series_2d = evecs(:, cmp)' * eeg_data_2d;
        ged_time_series_filtered_2d = evecs(:, cmp)' * eeg_data_filtered_2d;

        % Compute power for filtered signal
        ged_time_series_filtered_2d_power = abs(hilbert(ged_time_series_filtered_2d)) .^ 2;

        % Reshape time series to 3d
        ged_time_series(cmp, :, :) = reshape(ged_time_series_2d, [length(EEG.times), EEG.trials]);
        ged_time_series_filtered(cmp, :, :) = reshape(ged_time_series_filtered_2d_power, [length(EEG.times), EEG.trials]);

    end

    % Threshold eigenvalues
    suprathresh_eval_idx_p95 = find(evals > p95_threshold);

    % Threshold eigenvalues
    thresh_eval_sd = 1; % In sd
    thresh_eigenvalue = median(evals) + std(evals) * thresh_eval_sd;
    suprathresh_eval_idx = find(evals > thresh_eigenvalue);  

    % Determine electrode distances based on cartesian coordinates
    dists = zeros(EEG.nbchan);
    cart_coords = [cell2mat({EEG.chanlocs.X})', cell2mat({EEG.chanlocs.Y})', cell2mat({EEG.chanlocs.Z})'];
    for ch1 = 1 : EEG.nbchan
        crds1 = cart_coords(ch1, :);
        len1 = sqrt(sum(crds1 .^ 2));
        for ch2 = 1 : EEG.nbchan
            crds2 = cart_coords(ch2, :);
            len2 = sqrt(sum(crds2 .^ 2));
            if ch1 == ch2
                dists(ch1, ch2) = 0;
            else
                r = (len1 + len2) / 2; % Estimate sphere radius from len1 and len2
                theta = acos(dot(crds1, crds2) / (len1 * len2)); % Angle between A & B in radians
                dists(ch1, ch2) = r * theta; % Arc length = radius * theta
            end
        end
    end

    % Create spatial filter template
    focuschan = 65;
    template_topo = dists(focuschan, :) / max(dists(focuschan, :)); % Normalize distances
    template_topo = ones(size(template_topo)) - template_topo; % Invert

    % Find highest similarity in supra-threshold components
    cmpsim = 0;
    chosen_cmp = 0;
    for cmp = 1 : EEG.nbchan
        if ismember(cmp, suprathresh_eval_idx)
            tmp_cmp = ged_maps(cmp, :) / max(ged_maps(cmp, :)); % Normalize
            tmp = corrcoef(tmp_cmp, template_topo);
            tmp = abs(tmp(1, 2));
            if tmp > cmpsim
                cmpsim = tmp;
                chosen_cmp = cmp;
            end
        end	
    end

    if chosen_cmp ~= 0

        % Get selected component topography
        cmp_topo = ged_maps(chosen_cmp, :);

        % Save filter topography
        figure('Visible', 'off'); clf;
        topoplot(cmp_topo, EEG.chanlocs, 'electrodes', 'off', 'numcontour', 0)
        saveas(gcf, [PATH_OUT, 'selected_topo_' subject '.png']);
        
        % Get selected component signal
        cmp_time_series = squeeze(ged_time_series(chosen_cmp, :, :));
        cmp_time_series_filtered = squeeze(ged_time_series_filtered(chosen_cmp, :, :));

        % Create output struct
        ged_data = struct();
        ged_data.subject = subject;
        ged_data.ged_time_series = cmp_time_series;
        ged_data.ged_time_series_filtered = cmp_time_series_filtered;
        ged_data.ged_maps = ged_maps;
        ged_data.times = EEG.times;
        ged_data.trialinfo = EEG.trialinfo;

        % Save output struct
        save([PATH_OUT, subject, '_ged_data.mat'], 'ged_data')

        % Plot filter topographies
        figure('Visible', 'off');

        % Iterate first 5
        for cmp = 1 : 5

            title_str = '';

            % Plot filter topo
            subplot(5, 4, (cmp - 1) * 4 + 1)
            cmp_topo = ged_maps(cmp, :);
            topoplot(cmp_topo, EEG.chanlocs, 'electrodes', 'off', 'numcontour', 0);
            if cmp == chosen_cmp
                title_str = [title_str, 's'];
            end
            if evals(cmp) > p95_threshold
                title_str = [title_str, '*'];
            end
            if evals(cmp) > thresh_eigenvalue
                title_str = [title_str, '+'];
            end
            title(title_str)

            % Plot component average
            subplot(5, 4, (cmp - 1) * 4 + 2)
            pd = mean(squeeze(ged_time_series(cmp, :, :)), 2);  
            plot(EEG.times, pd, 'm', 'LineWidth', 1.5)

            % Plot component average theta
            subplot(5, 4, (cmp - 1) * 4 + 3)
            pd = mean(squeeze(ged_time_series_filtered(cmp, :, :)), 2);
            plot(EEG.times, pd, 'm', 'LineWidth', 1.5)

            % Plot eigenspectrum
            subplot(5, 4, (cmp - 1) * 4 + 4)
            plot([1:EEG.nbchan], evals, 'm*', 'MarkerSize', 6)
            yline(p95_threshold, 'k')
            yline(thresh_eigenvalue, 'g')
            xlim([0, EEG.nbchan + 1])
            title('eigenspectrum')

        end
        saveas(gcf, [PATH_OUT, 'spatial_filter_' subject '.png']);

    end

end





