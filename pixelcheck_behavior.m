clear all

% Paths 
PATH_EEGLAB = '/home/plkn/eeglab2025.0.0/';
PATH_RAW = '/mnt/data_dump/pixelcheck/0_eeg/';
PATH_BEHAVIOR = '/mnt/data_dump/pixelcheck/1_behavioral_data/';

% Init EEGlab
addpath(PATH_EEGLAB);
eeglab;

subject_list = {'VP001', 'VP003', 'VP004', 'VP005', 'VP006', 'VP007', 'VP009', 'VP010', 'VP011', ...)
                'VP012', 'VP013', 'VP014', 'VP015', 'VP016', 'VP017', 'VP018', 'VP019', 'VP020', 'VP021', ...
                'VP022', 'VP023', 'VP024', 'VP025', 'VP026', 'VP027', 'VP028', 'VP029', 'VP030', 'VP0301', ...
                'VP032', 'VP033_1', 'VP034', 'VP035', 'VP036', 'VP037', 'VP038', 'VP041'};

subject_list = {'VP004'};

% Init result array
behavioral_data = [];
behav_counter = 0;

% Loop subjects
for s = 1 : length(subject_list)

        % Get id
    if strcmp(subject_list{s}, 'VP0301')
        id = 31;
    else
        id = str2double(subject_list{s}(3 : 5));
    end

    if ismember(id, [35, 36, 37])

        % Load raw data
        EEG = pop_loadbv(PATH_RAW, ['Pixelcheck_', subject_list{s}, '.vhdr'], [], 1);

    else

        % Load raw data
        EEG = pop_loadbv(PATH_RAW, ['PixelCheck_', subject_list{s}, '.vhdr'], [], 1);
    
    end

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
    trialinfo = array2table(trialinfo, 'VariableNames', {'id','trial_nr','block_nr','reward_cue','reward_condition','cue_target_latency','correct_key_color','ma_condition','fb_correct','fb_flipped', 'last_fb_flipped','accuracy','rt'});
    aa=bb
    % Code trialnumber in block
    trialinfo.trial_in_block = mod(trialinfo.trial_nr, 120);
    trialinfo.trial_in_block(trialinfo.trial_in_block == 0) = 120;

    if s == 1
        T = trialinfo;
    else
        T = [T; trialinfo];
    end

    % Set trial exclusion criteria
    %idx_keep = EEG.trialinfo.trial_in_block > 40 & EEG.trialinfo.block_nr >=2;

    % Exclude trials
    %trialinfo = EEG.trialinfo(idx_keep, :);

    % Loop conditions
    % for agency = 1 : 3
    %     for reward = 0 : 1

    %         behav_counter = behav_counter + 1;

    %         % Get condition idx
    %         idx = trialinfo.ma_condition == agency & trialinfo.reward_condition == reward;

    %         % Get average rt for condition
    %         rt = mean(trialinfo.rt(trialinfo.accuracy == 1 & idx));

    %         % Get accuracy for condition
    %         accuracy = sum(trialinfo.accuracy == 1 & idx) / sum(idx);

    %         behavioral_data(behav_counter, :) = [trialinfo.id(1, 1), agency, reward, rt, accuracy];
    %     end
    % end

end % end subject loop

% Cast table
%behavioral_data = array2table(behavioral_data, 'VariableNames', {'id', 'agency', 'reward', 'rt', 'accuracy'});

% Save
writetable(T, [PATH_BEHAVIOR, 'behavioral_data.csv']);
