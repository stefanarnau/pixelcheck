clear all

% Paths 
PATH_EEGLAB = '/home/plkn/eeglab2025.0.0/';
PATH_RAW = '/mnt/data_dump/pixelcheck/0_eeg/';
PATH_AUTOCLEANED_ERP = '/mnt/data_dump/pixelcheck/2_cleaned_cuelocked_erp/';
PATH_AUTOCLEANED_TF = '/mnt/data_dump/pixelcheck/2_cleaned_cuelocked_tf/';

% Subject list 
subject_list = {'VP001', 'VP003', 'VP004', 'VP005', 'VP006', 'VP007', 'VP009', 'VP010', 'VP011', ...)
                'VP012', 'VP013', 'VP014', 'VP015', 'VP016', 'VP017', 'VP018', 'VP019', 'VP020', 'VP021', ...
                'VP022', 'VP023', 'VP024', 'VP025', 'VP026', 'VP027', 'VP028', 'VP029', 'VP030', 'VP0301', ...
                'VP032', 'VP033_1', 'VP034', 'VP035', 'VP036', 'VP037', 'VP038', 'VP041'};

% Initialize EEGLab
addpath(PATH_EEGLAB);
eeglab;
channel_location_file = which('standard-10-5-cap385.elp');

% Iterate subjects
for s = 1 : length(subject_list)

    % Get id
    if strcmp(subject_list{s}, 'VP0301')
        id = 31;
    else
        id = str2double(subject_list{s}(3 : 5));
    end

    if ismember(id, [35, 36, 37])

        % Load raw data
        EEG = pop_loadbv(PATH_RAW, ['Pixelcheck_', subject_list{s}, '.vhdr'], [], []);

    else

        % Load raw data
        EEG = pop_loadbv(PATH_RAW, ['PixelCheck_', subject_list{s}, '.vhdr'], [], []);
    
    end

    % STEP 1: Event coding ==================================================================================================

    trialinfo = [];
    counter = 0;
    last_fb_flipped = 0;
    last_block = 1;
    for e = 1 : length(EEG.event)
        if strcmpi(EEG.event(e).type(1), 'S')

            % Get event number from 'S' triggers
            enum = str2double(EEG.event(e).type(2:end));

            % If event number in range, then this is grid onset
            if enum >= 7 && enum <= 78

                % Increase counter...
                counter = counter + 1;
                
                % Loop backwards to find corresponding cue onset
                f = e;
                cuenum = 0;
                while ~ismember(cuenum, [91, 92])
                    f = f - 1;
                    cuenum = str2double(EEG.event(f).type(2:end));
                end

                % Mark cue onsets
                EEG.event(f).type = 'cue';
                EEG.event(f).code = 'cue';

                % Get reward condition
                reward_cue = (cuenum == 91) + 2*(cuenum == 92); % X=1, O=2
                if mod(id, 2) == 0 % If id even
                    reward_condition = reward_cue - 1; % If X->0, if O->1
                else
                    reward_condition = 2 - reward_cue; % If X->1, if O->0
                end

                % Get cue target interval length in frames
                cue_target_latency = EEG.event(e).latency - EEG.event(f).latency;

                % Get correct key color
                correct_key_color = mod(enum, 2) == 0 + 1;

                % Get block number
                block_nr = ceil(enum / 6) - 1;

                % Check if block changed
                if block_nr ~= last_block
                    last_fb_flipped = 0;
                    last_block = block_nr;
                end

                % Get ma condition
                cond_map = [1, 1, 2, 2, 3, 3];
                ma_condition = cond_map(mod(enum - 1, 6) + 1);

                % Loop forward for feedback
                f = e;
                fbnum = 0;
                while ~ismember(fbnum, 81:89) & f < length(EEG.event)
                    f = f + 1;
                    fbnum = str2double(EEG.event(f).type(2:end));
                end

                % Mark feedback events
                EEG.event(f).type = 'feedback';
                EEG.event(f).code = 'feedback';

                % Get feedback type
                if ismember(fbnum, [81, 83, 86])
                    acc = 1; fb_correct = 1; fb_flipped = 0;
                elseif ismember(fbnum, [82, 84, 87])
                    acc = 0; fb_correct = 0; fb_flipped = 0;
                elseif ismember(fbnum, [85, 88])
                    acc = 1; fb_correct = 0; fb_flipped = 1;
                else
                    acc = 0; fb_correct = 2; fb_flipped = 0; 
                end

                % Get response times
                rt = EEG.event(f).latency - EEG.event(e).latency;

                % Fill table
                trialinfo(counter,:) = [id,...
                                        counter,...
                                        block_nr,...
                                        reward_cue,...
                                        reward_condition,...
                                        cue_target_latency,...
                                        correct_key_color,...
                                        ma_condition,...
                                        fb_correct,...
                                        fb_flipped,...
                                        last_fb_flipped,...
                                        acc,...
                                        rt];

                last_fb_flipped = fb_flipped;
            end
        end
    end
        
    % Make table
    EEG.trialinfo = array2table(trialinfo, 'VariableNames', {'id','trial_nr','block_nr','reward_cue','reward_condition','cue_target_latency','correct_key_color','ma_condition','fb_correct','fb_flipped', 'last_fb_flipped','accuracy','rt'});

    % STEP 2: Channel structure ==================================================================================================

    % Remove sync channel
    EEG = pop_select(EEG, 'nochannel', {'sync'});

    % Add FCz (no here, as it has been used as online reference)
    EEG.data(end + 1, :) = 0;
    EEG.nbchan = size(EEG.data, 1);
    EEG.chanlocs(end + 1).labels = 'FCz';
    EEG = pop_chanedit(EEG, 'lookup', channel_location_file);

    % Save channel locations
    EEG.chanlocs_original = EEG.chanlocs;

    % Rereference to CPz (Fill FCz with values, get rid of CPz, for now)
    EEG = pop_reref(EEG, 'CPz');

    % Channel rejection based on kurtosos and spectral criteria
    [EEG, rej1]    = pop_rejchan(EEG,    'elec', 1:EEG.nbchan, 'threshold', 5, 'norm', 'on', 'measure', 'kurt');
    [EEG, rej2]    = pop_rejchan(EEG,    'elec', 1:EEG.nbchan, 'threshold', 3, 'norm', 'on', 'measure', 'spec');
    EEG.chans_rejected_combined    = unique([rej1, rej2]);

    % Interpolate rejected channels (and CPz)
    EEG = pop_interp(EEG, EEG.chanlocs_original, 'spherical');

    % Rereference all channels to Common average reference
    EEG = pop_reref(EEG, []);

    % Resampling data for ICA and TF analyses
    ERP = EEG;
    EEG = pop_resample(EEG, 200);

    % Filter the data
    ERP = pop_basicfilter(ERP, 1:ERP.nbchan, 'Cutoff', [0.1, 30],    'Design', 'butter', 'Filter', 'bandpass', 'Order', 4);
    EEG = pop_basicfilter(EEG, 1:EEG.nbchan, 'Cutoff', [2, 30],    'Design', 'butter', 'Filter', 'bandpass', 'Order', 4);

    % Create epochs
    [ERP, epoch_idx_erp] = pop_epoch(ERP, {'cue'}, [-0.2, 2.5]);
    [EEG, epoch_idx_tf]  = pop_epoch(EEG, {'cue'}, [-0.8, 3]);

    % Reduce trialinfo
    ERP.trialinfo = ERP.trialinfo(epoch_idx_erp, :);
    EEG.trialinfo = EEG.trialinfo(epoch_idx_tf, :);

    % Remove baseline
    ERP = pop_rmbase(ERP, [-200, 0]);
    EEG = pop_rmbase(EEG, [-200, 0]);

    % Detect and reject bad epochs
    [EEG, EEG.rejected_epochs] = pop_autorej(EEG, 'nogui', 'on');

    % Remove those epochs also from trialinfo
    EEG.trialinfo(EEG.rejected_epochs, :) = [];

    % Run ICA
    EEG = pop_runica(EEG, 'extended', 1, 'interrupt', 'on', 'PCA', 64 - (length(EEG.chans_rejected_combined) + 2));

    % Run IC Label
    EEG = iclabel(EEG);

    % Define non brain ICs
    EEG.nobrainer = find(EEG.etc.ic_classification.ICLabel.classifications(:, 1) < 0.3 | EEG.etc.ic_classification.ICLabel.classifications(:, 3) > 0.3);

    % Copy ICs to erpset
    ERP = pop_editset(ERP, 'icachansind', 'EEG.icachansind', 'icaweights', 'EEG.icaweights', 'icasphere', 'EEG.icasphere');
    ERP.etc = EEG.etc;
    ERP.nobrainer = EEG.nobrainer;

    % Remove non-brain components
    ERP = pop_subcomp(ERP, ERP.nobrainer, 0);
    EEG = pop_subcomp(EEG, EEG.nobrainer, 0);

    % Detect and reject bad epochs in erp set
    [ERP, ERP.rejected_epochs] = pop_autorej(ERP, 'nogui', 'on');

    % Remove those epochs also from trialinfo
    ERP.trialinfo(ERP.rejected_epochs, :) = [];

    % STEP 6: Save data ==================================================================================================

    % Save erp set
    pop_saveset(ERP, 'filename', [num2str(id), '_cleaned.set'], 'filepath', PATH_AUTOCLEANED_ERP);

    % Save tf set
    pop_saveset(EEG, 'filename', [num2str(id), '_cleaned.set'], 'filepath', PATH_AUTOCLEANED_TF);

end % End subject iteration
