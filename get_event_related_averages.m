function[erp, erp_time, ersp_raw, ersp_db_general, ersp_db_specific, itpc, tf_time, tf_frqs] = get_event_related_averages(EEG, fserp, fstf, cndidx, varargin)

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
	% t       	       : A vector containing the latencies of a trial (e.g. EEG.times in an epoched dataset)
	% fs     	       : The sampling rate (e.g. EEG.srate)
    % cndidx           : Trial indices for each condition
    %
	% optional: 
	% 'blerp         ' : Latencies to use for the erp baseline. Default is [-200 0].
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
	% erp              : The erp matrix as channel x erp
	% erp_time         : The time vector for erp data.
    % ersp_raw         : The ersp matrix in mV^2 values as channel x frequency x time
	% ersp_db_general  : The ersp matrix in decibel values using a condition-general baseline as channel x frequency x time
    % ersp_db_specific : The ersp matrix in decibel values using a condition-specific baseline as channel x frequency x time
	% itpc             : The itpc matrix as channel x frequency x time
	% tf_time          : The time vector for time-frequency data. Differs from erp_time as edge artifacts are cropped.
	% tf_frqs          : The frequency vector for time-frequency data.
	%
	% ------------------------
	% Stefan Arnau, 27.03.2019
	% ------------------------

	% Check if any input args
	if nargin < 3
        error('Not enough input arguments... :)');
        return;
	end

	% Init input parser
	p = inputParser;

	% Set Defaults
	default_blerp = [-200, 0];
	default_blersp = [-500, -200];
	default_n_frq = 20;
	default_frqrange = [2, 30];
	default_fwhmtrange = [500, 100];

	% Parse inputs and set defaults
	p.FunctionName  = mfilename;
	p.CaseSensitive = false;
	p.addRequired('EEG', @isstruct);
	p.addRequired('fserp', @isnumeric);
    p.addRequired('fstf', @isnumeric);
    p.addRequired('cndidx', @iscell);
	p.addParamValue('blerp', default_blerp, @isnumeric);
	p.addParamValue('blersp', default_blersp, @isnumeric);
	p.addParamValue('n_frq', default_n_frq, @isnumeric);
	p.addParamValue('frqrange', default_frqrange, @isnumeric);
	p.addParamValue('fwhmtrange', default_fwhmtrange, @isnumeric);
    parse(p, EEG, fserp, fstf, cndidx, varargin{:});
    
    % Resample
    DERP = pop_resample(p.Results.EEG, p.Results.fserp);
    DTF = pop_resample(p.Results.EEG, p.Results.fstf);
    pruneidx = dsearchn(DTF.times', [DTF.times(1) + 300, DTF.times(end) - 300]');
    tf_time = DTF.times(pruneidx(1) : pruneidx(2));
    erp_time = DERP.times;

    % Resmats
    erp = zeros(length(p.Results.cndidx), size(DERP.data, 1), size(DERP.data, 2));
    ersp = zeros(length(p.Results.cndidx), size(DTF, 1), p.Results.n_frq, length(tf_time));
    itpc = zeros(length(p.Results.cndidx), size(DTF, 1), p.Results.n_frq, length(tf_time));

    % Iterate conditions
    for cnd = 1 : length(p.Results.cndidx)
        
        % Calc erps
        d = DERP.data(:, :, p.Results.cndidx{cnd});
        blerpidx = dsearchn(DERP.times', p.Results.blerp');
        for c = 1 : size(d, 1)
            fprintf('\ncondition %i/%i | Calculating ERP chan %i/%i...', cnd, length(p.Results.cndidx), c, size(d, 1));
            erp(cnd, c, :) = mean(bsxfun(@minus, squeeze(d(c, :, :))', mean(squeeze(d(c, blerpidx(1) : blerpidx(2), :))', 2)), 1);
            fprintf('  done');
        end

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

    % Talk more
    fprintf('  done\n\n');

end


