function[ersp_raw, ersp_db_general, ersp_db_specific, itpc, tf_time, tf_frqs] = get_event_related_tf(EEG, cndidx, varargin)

	%
	%
	% WHAT THIS FUNCTION DOES:
	% ------------------------
	%
    % This function takes eeg data and returns erp, db-baseline normalized ersp and itpc for each channel, as well 
    % as the corresponding time and frequency vectors used for time-frequency decomposition.
	%
	% USAGE: [erp, erp_time, ersp, itpc, tf_time, tf_frqs] = get_tf_stuff_from_EEG(EEG, 500, 200, {[1,2,3,4,5], [6,7,8,9,17]});
	%
	%
	% INPUT ARGUMENTS:
	% ----------------
	%
	% needed: 
	% EEG              : The EEG struct
    % cndidx           : Trial indices for each condition
    %
	% optional: 
	% 'blersp'         : Latencies to use for the ersp baseline. Default is [-500 -200].
	% 'n_frq     '     : Number of frequencies to use in tf-decomposition. Default is 20.
	% 'frqrange    '   : Frequency range for tf-decomposition. Default is [3, 30].
    % 'fwhmrrange'     : Range of 'full-width-at-half-maximum' in the time domain (i.e. ms-values) 
    %                    for wavelet construction. Default is [400, 100].
	%
	%
	% OUTPUT ARGUMENTS:
	% -----------------
	%	
    % ersp_raw         : The ersp matrix in mV^2 values as channel x frequency x time
	% ersp_db_general  : The ersp matrix in decibel values using a condition-general baseline as channel x frequency x time
    % ersp_db_specific : The ersp matrix in decibel values using a condition-specific baseline as channel x frequency x time
	% itpc             : The itpc matrix as channel x frequency x time
	% tf_time          : The time vector for time-frequency data. Differs from erp_time as edge artifacts are cropped.
	% tf_frqs          : The frequency vector for time-frequency data.
	%
	% ------------------------
	% Stefan Arnau, 25.11.2025
	% ------------------------

	% Init input parser
	p = inputParser;

	% Set Defaults
	default_blersp = [-500, -200];
	default_n_frq = 20;
	default_frqrange = [2, 30];
	default_fwhmtrange = [500, 100];

	% Parse inputs and set defaults
	p.FunctionName  = mfilename;
	p.CaseSensitive = false;
	p.addRequired('EEG', @isstruct);
    p.addRequired('cndidx', @iscell);
	p.addParamValue('blersp', default_blersp, @isnumeric);
	p.addParamValue('n_frq', default_n_frq, @isnumeric);
	p.addParamValue('frqrange', default_frqrange, @isnumeric);
	p.addParamValue('fwhmtrange', default_fwhmtrange, @isnumeric);
    parse(p, EEG, cndidx, varargin{:});
    
    % Set prune time
    DTF = EEG;
    pruneidx = dsearchn(DTF.times', [DTF.times(1) + 500, DTF.times(end) - 500]');
    tf_time = DTF.times(pruneidx(1) : pruneidx(2));

    % Resmats
    ersp = zeros(length(p.Results.cndidx), EEG.nbchan, p.Results.n_frq, length(tf_time));
    itpc = zeros(length(p.Results.cndidx), EEG.nbchan, p.Results.n_frq, length(tf_time));

    % Iterate conditions
    for cnd = 1 : length(p.Results.cndidx)

        % If too few trials
        if length(p.Results.cndidx{cnd}) < 3
                ersp(cnd, :, :, :) = nan(EEG.nbchan, p.Results.n_frq, length(tf_time));
                itpc(cnd, :, :, :) = nan(EEG.nbchan, p.Results.n_frq, length(tf_time));
        else

            % Set some variables
            d = DTF.data(:, :, p.Results.cndidx{cnd});
            wtime = -2 : 1 / DTF.srate : 2;
            halfw = (length(wtime) - 1) / 2;
            nconv = size(d, 2) * size(d, 3) + length(wtime) - 1;
            tf_frqs = logspace(log10(p.Results.frqrange(1)), log10(p.Results.frqrange(2)), p.Results.n_frq);
            fwhms = logspace(log10(p.Results.fwhmtrange(1)), log10(p.Results.fwhmtrange(2)), p.Results.n_frq);
            cmwX = zeros(p.Results.n_frq, nconv);

            % Build wavelets
            fprintf('\ncondition %i/%i | Building wavelets for tf-analysis...', cnd, length(p.Results.cndidx));
            for frq = 1 : p.Results.n_frq
                cmw = exp(2 * 1i * pi * tf_frqs(frq) .* wtime) .* exp((-4 * log(2) * wtime.^2) ./ (fwhms(frq) / 1000)^2); % Build
                cmw = fft(cmw ./ max(cmw), nconv); % Normalize in time domain and fft
                cmwX(frq, :) = cmw ./ max(cmw); % Normalize in frq domain
            end
            fprintf('  done');

            % Get tf trial mesures
            for c = 1 : size(d, 1)
                fprintf('\ncondition %i/%i | Calculating time-frequency decomposition chan %i/%i...', cnd, length(p.Results.cndidx), c, size(p.Results.EEG.data, 1));

                % Convolute
                pow = zeros(p.Results.n_frq, size(d, 2), size(d, 3));
                pha = zeros(p.Results.n_frq, size(d, 2), size(d, 3)); 
                for frq = 1 : p.Results.n_frq
                    dX = fft(reshape(double(squeeze(d(c, :, :))), 1, []), nconv);
                    as = ifft(cmwX(frq, :) .* dX);
                    as = as(halfw + 1 : end - halfw);
                    as = reshape(as, size(d, 2), size(d, 3));
                    pow(frq, :, :) = abs(as) .^ 2;
                    pha(frq, :, :) = angle(as);
                end

                % Cut edge artifacts
                pow = pow(:, pruneidx(1) : pruneidx(2), :);
                pha = pha(:, pruneidx(1) : pruneidx(2), :);

                % Get tf measures
                ersp(cnd, c, :, :) = squeeze(mean(pow, 3));
                itpc(cnd, c, :, :) = abs(squeeze(mean(exp(1i*pha), 3)));

                fprintf('  done');
            end
        end
    end

    % New matrices are needed
    ersp_raw = ersp;
    ersp_db_general = ersp;
    ersp_db_specific = ersp;

    % Apply condition general baseline
    fprintf('\nBaselining data...');
    blerspidx = dsearchn(tf_time', p.Results.blersp');
    for c = 1 : DTF.nbchan

        % Get condition-general baseline values
        blvals_general = mean(squeeze(mean(squeeze(ersp(:, c, :, blerspidx)), 1)), 2);

        % Iterate conditions
        for cnd = 1 : length(p.Results.cndidx)

            % If too few trials
            if length(p.Results.cndidx{cnd}) < 3
                    ersp_db_general(cnd, c, :, :) = nan(p.Results.n_frq, length(tf_time));
                    ersp_db_general(cnd, c, :, :) = nan(p.Results.n_frq, length(tf_time));
            else

                % Get condition-specific baseline values
                blvals_specific = mean(squeeze(ersp(cnd, c, :, blerspidx)), 2);

                % Get condition x channel data
                tmp = squeeze(ersp(cnd, c, :, :));

                % Apply condition-general baseline
                ersp_db_general(cnd, c, :, :) = 10 * log10(bsxfun(@rdivide, tmp, blvals_general));

                % Apply condition-specific baseline
                ersp_db_specific(cnd, c, :, :) = 10 * log10(bsxfun(@rdivide, tmp, blvals_specific));

            end

        end
    end

    % Talk more
    fprintf('  done\n\n');

end


