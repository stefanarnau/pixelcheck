%% ======= Time frequency analysis for feedback locked data ===========
% ========= 15.07.2025 =================
% Based off pixelcheck_tf_cluster skript=======


clear all;


% PATH FOR MAC
  PATH_EEGLAB =  '/Users/grote/Library/CloudStorage/OneDrive-ifado.de/Dokumente/eeglab2025.0.0/';
  PATH_FIELDTRIP = '/Users/grote/Library/CloudStorage/OneDrive-ifado.de/Dokumente/MATLAB/fieldtrip';
  PATH_AUTOCLEANED = '/Users/grote/Library/CloudStorage/OneDrive-ifado.de/Dokumente/MATLAB/Pixelcheck/feedback_locked/';
  PATH_TF_DATA = '/Users/grote/Library/CloudStorage/OneDrive-ifado.de/Dokumente/MATLAB/Pixelcheck/tf_feedback/';
  PATH_FIRST_PLOTS = '/Users/grote/Library/CloudStorage/OneDrive-ifado.de/Dokumente/MATLAB/Pixelcheck/cluster_plots_tf/feedback/';


% Subject list 
subject_list = {'VP001', 'VP003', 'VP004', 'VP005', 'VP006', 'VP007', 'VP009', 'VP010', 'VP011', ...)
                'VP012', 'VP013', 'VP014', 'VP015', 'VP016', 'VP017', 'VP018', 'VP019', 'VP020', 'VP021', ...
                'VP022', 'VP023', 'VP024', 'VP025', 'VP026', 'VP027', 'VP028', 'VP029', 'VP030', 'VP031', ...
                'VP032', 'VP033', 'VP034', 'VP035', 'VP036', 'VP037', 'VP038', 'VP041'}; % 

% Init eeglab
addpath(PATH_EEGLAB);
eeglab;

% Init fieldtrip
addpath(PATH_FIELDTRIP);
ft_defaults;

% SWITCH: Switch parts of script on/off
to_execute = {'part1'}; %toggle to part 2 , 'part2'

% Part 1: Calculate ersp & itpc for 12 ERP conditions
if ismember('part1', to_execute)

    % Lade Info von erstem Subjekt
    EEG = pop_loadset('filename', [subject_list{1} '_cleaned_feedback_tf_.set.set'], 'filepath', PATH_AUTOCLEANED, 'loadmode', 'info');

    % Wavelet Parameter
    n_frq = 30;
    frqrange = [2, 30];
    tfres_range = [600, 240];
    wtime = -2 : 1/EEG.srate : 2;

    hz = linspace(0, EEG.srate, length(wtime));

    tf_freqs = logspace(log10(frqrange(1)), log10(frqrange(2)), n_frq);
    fwhmTs = logspace(log10(tfres_range(1)), log10(tfres_range(2)), n_frq);

    % Erzeuge Wavelets
    cmw = zeros(length(tf_freqs), length(wtime));
    for frq = 1:length(tf_freqs)
        cmw(frq,:) = exp(2*1i*pi*tf_freqs(frq).*wtime) .* exp((-4*log(2)*wtime.^2) ./ (fwhmTs(frq)/1000)^2);
        cmw(frq,:) = cmw(frq,:) ./ max(cmw(frq,:));
    end

    % Analysefenster
    prune_times = [-500, 2500];
    tf_times = EEG.times(dsearchn(EEG.times', prune_times(1)):dsearchn(EEG.times', prune_times(2)));

    % Loop subjects
    for s = 1:length(subject_list)

        subject = subject_list{s};
        fprintf('\n== Subjekt %s ==\n', subject);

        EEG = pop_loadset('filename', [subject '_cleaned_feedback_tf_.set.set'], 'filepath', PATH_AUTOCLEANED);

        % Trials aus Blöcken 2–12
        EEG.trialinfo.trial_in_block = mod(EEG.trialinfo.trial_nr,120);
        EEG.trialinfo.trial_in_block(EEG.trialinfo.trial_in_block==0)=120;

        idx_keep = EEG.trialinfo.trial_in_block > 40 & EEG.trialinfo.block_nr >=2;
        EEG.data = EEG.data(:,:,idx_keep);
        EEG.trialinfo = EEG.trialinfo(idx_keep,:);
        EEG.trials = sum(idx_keep);

        eeg_data = double(EEG.data);
        chanlocs = EEG.chanlocs;

        nchan = EEG.nbchan;
        ntime = length(tf_times);
        nfreq = length(tf_freqs);

        ersps = single(zeros(12, nchan, nfreq, ntime));
        itpcs = single(zeros(12, nchan, nfreq, ntime));

        % Schleife über Kanäle
        for ch = 1:nchan

            fprintf('  Kanal %d/%d\n', ch, nchan);

            signal = squeeze(eeg_data(ch,:,:));

            convlen = size(signal,1) * size(signal,2) + size(cmw,2)-1;
            tmp = fft(reshape(signal,1,[]), convlen);

            powcube = NaN(nfreq, EEG.pnts, EEG.trials);
            phacube = NaN(nfreq, EEG.pnts, EEG.trials);

            cmwX = zeros(nfreq, convlen);
            for f = 1:nfreq
                cmwX(f,:) = fft(cmw(f,:), convlen);
                cmwX(f,:) = cmwX(f,:) ./ max(cmwX(f,:));
            end

            for f = 1:nfreq
                as = ifft(cmwX(f,:) .* tmp);
                as = as(((size(cmw,2)-1)/2)+1:end-((size(cmw,2)-1)/2));
                as = reshape(as, EEG.pnts, EEG.trials);
                powcube(f,:,:) = abs(as).^2;
                phacube(f,:,:) = exp(1i*angle(as));
            end

            % Crop Zeitfenster
            idx_t1 = dsearchn(EEG.times', prune_times(1));
            idx_t2 = dsearchn(EEG.times', prune_times(2));

            powcube = powcube(:,idx_t1:idx_t2,:);
            phacube = phacube(:,idx_t1:idx_t2,:);

            % Baseline [-500:-200]
            [~, bl1] = min(abs(tf_times - (-500)));
            [~, bl2] = min(abs(tf_times - (-200)));
            blvals = squeeze(mean(mean(powcube(:,bl1:bl2,:),3),2));

            % Schleife über Bedingungen 1–12
            for cond = 1:12
                switch cond
                    % SELF
                    case 1
                        idx = EEG.trialinfo.ma_condition==2 & EEG.trialinfo.fb_correct==0 & EEG.trialinfo.fb_flipped==0 & EEG.trialinfo.reward_condition==0;
                    case 2
                        idx = EEG.trialinfo.ma_condition==2 & EEG.trialinfo.fb_correct==0 & EEG.trialinfo.fb_flipped==0 & EEG.trialinfo.reward_condition==1;
                    case 3
                        idx = EEG.trialinfo.ma_condition==2 & EEG.trialinfo.fb_correct==0 & EEG.trialinfo.fb_flipped==1 & EEG.trialinfo.reward_condition==0;
                    case 4
                        idx = EEG.trialinfo.ma_condition==2 & EEG.trialinfo.fb_correct==0 & EEG.trialinfo.fb_flipped==1 & EEG.trialinfo.reward_condition==1;
                    case 5
                        idx = EEG.trialinfo.ma_condition==2 & EEG.trialinfo.fb_correct==1 & EEG.trialinfo.reward_condition==0;
                    case 6
                        idx = EEG.trialinfo.ma_condition==2 & EEG.trialinfo.fb_correct==1 & EEG.trialinfo.reward_condition==1;
                    % OTHER
                    case 7
                        idx = EEG.trialinfo.ma_condition==3 & EEG.trialinfo.fb_correct==0 & EEG.trialinfo.fb_flipped==0 & EEG.trialinfo.reward_condition==0;
                    case 8
                        idx = EEG.trialinfo.ma_condition==3 & EEG.trialinfo.fb_correct==0 & EEG.trialinfo.fb_flipped==0 & EEG.trialinfo.reward_condition==1;
                    case 9
                        idx = EEG.trialinfo.ma_condition==3 & EEG.trialinfo.fb_correct==0 & EEG.trialinfo.fb_flipped==1 & EEG.trialinfo.reward_condition==0;
                    case 10
                        idx = EEG.trialinfo.ma_condition==3 & EEG.trialinfo.fb_correct==0 & EEG.trialinfo.fb_flipped==1 & EEG.trialinfo.reward_condition==1;
                    case 11
                        idx = EEG.trialinfo.ma_condition==3 & EEG.trialinfo.fb_correct==1 & EEG.trialinfo.reward_condition==0;
                    case 12
                        idx = EEG.trialinfo.ma_condition==3 & EEG.trialinfo.fb_correct==1 & EEG.trialinfo.reward_condition==1;
                end

                pow_this = mean(powcube(:,:,idx),3);
                pha_this = mean(phacube(:,:,idx),3);

                ersps(cond,ch,:,:) = 10*log10(bsxfun(@rdivide, pow_this, blvals));
                itpcs(cond,ch,:,:) = abs(pha_this);
            end
        end % Kanal

        % Speichern
        save([PATH_TF_DATA, subject '_ersps_feedback_12cond.mat'], 'ersps', 'tf_freqs', 'tf_times', 'chanlocs', '-v7.3');
        save([PATH_TF_DATA, subject '_itpcs_feedback_12cond.mat'], 'itpcs', 'tf_freqs', 'tf_times', 'chanlocs', '-v7.3');

    end % Subjekt

end % Part1


%%  PART 2 für ANOVA 

% clearvars;
% 
% % === Pfade & Settings ===
% PATH_TF_DATA = '/Users/grote/Library/CloudStorage/OneDrive-ifado.de/Dokumente/MATLAB/Pixelcheck/tf_feedback/';
% subject_list = {'VP001', 'VP003', 'VP004', 'VP005', 'VP006', 'VP007', 'VP009', 'VP010', 'VP011', ...
%                 'VP012', 'VP013', 'VP014', 'VP015', 'VP016', 'VP017', 'VP018', 'VP019', 'VP020', 'VP021', ...
%                 'VP022', 'VP023', 'VP024', 'VP025', 'VP026', 'VP027', 'VP028', 'VP029', 'VP030', 'VP031', ...
%                 'VP032', 'VP033', 'VP034', 'VP035', 'VP036', 'VP037', 'VP038', 'VP041'};
% 
% n_cond = 12;
% 
% % ROI definieren
% freq_roi = [4 7];                 % Theta, z.B. 4–7 Hz
% time_roi = [200 400];            % Zeitfenster in ms
% chan_roi = {'Fz', 'FCz', 'Cz'};  % Frontozentral
% 
% % === Tabelle vorbereiten ===
% tbl = table();
% tbl.Subject = zeros(length(subject_list),1);
% 
% % === Schleife über alle Probanden ===
% for s = 1:length(subject_list)
% 
%     subject = subject_list{s};
%     fprintf('\n📄 Lade %s …\n', subject);
% 
%     load(fullfile(PATH_TF_DATA, [subject '_ersps_feedback_12cond.mat'])); % lädt: ersps, tf_freqs, tf_times, chanlocs
% 
%     % Finde ROI-Indizes beim ersten Probanden
%     if s == 1
%         freq_idx = find(tf_freqs >= freq_roi(1) & tf_freqs <= freq_roi(2));
%         time_idx = find(tf_times >= time_roi(1) & tf_times <= time_roi(2));
%         chan_idx = find(ismember({chanlocs.labels}, chan_roi));
%     end
% 
%     tbl.Subject(s) = str2double(subject(3:5));
% 
%     % Schleife über Bedingungen
%     for c = 1:n_cond
%         roi_vals = mean(ersps(c,chan_idx,freq_idx,time_idx), [2 3 4], 'omitnan');
%         colname = sprintf('Cond%d',c);
%         tbl.(colname)(s,1) = roi_vals;
%     end
% end
% 
% % === Speichern ===
% outfile = fullfile(PATH_TF_DATA,'TF_ROI_JASP_Wide_12conditions.csv');
% writetable(tbl, outfile);
% 
% fprintf('ERSP ROI-Tabelle für JASP mit 12 Bedingungen gespeichert:\n%s\n', outfile);
% disp(tbl);

% %%  PART 2b: ITPC ROI-Tabelle für ANOVA (JASP-ready)
% 
% % Pfade und Settings
% PATH_TF_DATA = '/Users/grote/Library/CloudStorage/OneDrive-ifado.de/Dokumente/MATLAB/Pixelcheck/tf_feedback/';
% n_cond = 12;
% 
% % ROI definieren
% freq_roi = [4 7];         % Hz
% time_roi = [200 400];     % ms
% chan_roi = {'Fz', 'FCz', 'Cz'};  % Frontozentral
% 
% % Lade Meta-Daten aus der ersten Datei
% load(fullfile(PATH_TF_DATA, sprintf('%s_itpcs_feedback_12cond.mat', subject_list{1})));  % lädt itpcs, tf_freqs, tf_times, chanlocs
% 
% n_subj = length(subject_list);
% 
% freq_idx = find(tf_freqs >= freq_roi(1) & tf_freqs <= freq_roi(2));
% time_idx = find(tf_times >= time_roi(1) & tf_times <= time_roi(2));
% chan_idx = find(ismember({chanlocs.labels}, chan_roi));
% 
% roi_mean = @(data) squeeze(mean(data(chan_idx, freq_idx, time_idx), [1 2 3]));  % Mittelwert über ROI
% 
% % Tabelle initialisieren
% tbl_itpc = table((1:n_subj)', 'VariableNames', {'Subject'});
% 
% % Schleife über alle 12 Bedingungen
% for c = 1:n_cond
%     vals = nan(n_subj,1);
% 
%     for s = 1:n_subj
%         subject = subject_list{s};
%         fname = fullfile(PATH_TF_DATA, sprintf('%s_itpcs_feedback_12cond.mat', subject));
%         load(fname);  % lädt itpcs, tf_freqs, tf_times, chanlocs
% 
%         % itpcs: [12 x nchan x nfreq x ntime]
%         roi_val = roi_mean(squeeze(itpcs(c,:,:,:)));
%         vals(s) = roi_val;
%     end
% 
%     colname = sprintf('Cond%d',c);
%     tbl_itpc.(colname) = vals;
% end
% 
% % Speichern
% writetable(tbl_itpc, fullfile(PATH_TF_DATA, 'TF_ITPC_ROI_JASP_Wide_12conditions.csv'));
% 
% disp('ITPC ROI-Tabelle für JASP mit 12 Bedingungen gespeichert:');
% disp(tbl_itpc);
 %% PLOTS
 %% ======================
% Einstellungen & Laden
% ======================

PATH_TF_DATA = '/Users/grote/Library/CloudStorage/OneDrive-ifado.de/Dokumente/MATLAB/Pixelcheck/tf_feedback/';

subject_list = {'VP001', 'VP003', 'VP004', 'VP005', 'VP006', 'VP007', 'VP009', 'VP010', 'VP011', ...
                'VP012', 'VP013', 'VP014', 'VP015', 'VP016', 'VP017', 'VP018', 'VP019', 'VP020', 'VP021', ...
                'VP022', 'VP023', 'VP024', 'VP025', 'VP026', 'VP027', 'VP028', 'VP029', 'VP030', 'VP031', ...
                'VP032', 'VP033', 'VP034', 'VP035', 'VP036', 'VP037', 'VP038', 'VP041'};

titles = {'Falsch Low', 'Falsch High', ...
          'Flipped Low', 'Flipped High', ...
          'Korrekt Low', 'Korrekt High'};

n_cond = 12;  % SELF (1–6) + OTHER (7–12)

mean_ersp = cell(1, n_cond);

% ======================
% Mittelwerte für alle 12 Bedingungen

for cond = 1:n_cond
    ersp_all = [];  % [nSubj x freq x time]
    for s = 1:length(subject_list)
        subj = subject_list{s};
        load(fullfile(PATH_TF_DATA, sprintf('%s_ersps_feedback_12cond.mat',subj)));
        
        chan_idx = find(strcmp({chanlocs.labels}, 'FCz'));  % ROI: FCz
        tmp = squeeze(ersps(cond,chan_idx,:,:));  % [freq x time]
        ersp_all(s,:,:) = tmp;

        if s==1
            tf_freqs = tf_freqs; tf_times = tf_times;  % speichern
        end
    end
    mean_ersp{cond} = squeeze(mean(ersp_all,1));  % [freq x time]
end

% ======================
% SELF Conditions (1–6) – TF Heatmaps bis 10Hz

figure;
for i = 1:6
    subplot(3,2,i);
    freq_mask = tf_freqs <= 10;  % nur bis 10 Hz
    data = mean_ersp{i}(freq_mask,:);
    contourf(tf_times, tf_freqs(freq_mask), data, 50, 'LineStyle', 'none');
    clim([-2, 2])
    colormap('jet')
    xlabel('Time (ms)'); ylabel('Freq (Hz)');
    title(titles{i});
    colorbar;
end
sgtitle('Time-Frequency (<=10Hz) – SELF conditions');

% ======================
% OTHER Conditions (7–12) – TF Heatmaps bis 10Hz

figure;
for i = 7:12
    subplot(3,2,i-6);
    freq_mask = tf_freqs <= 10;
    data = mean_ersp{i}(freq_mask,:);
    contourf(tf_times, tf_freqs(freq_mask), data, 50, 'LineStyle', 'none');
    clim([-2, 2])
    colormap('jet')
    xlabel('Time (ms)'); ylabel('Freq (Hz)');
    title(titles{i-6});
    colorbar;
end
sgtitle('Time-Frequency (<=10Hz) – OTHER conditions');

% ======================
% SELF Conditions (1–6) – TF Heatmaps bis 30Hz

figure;
for i = 1:6
    subplot(3,2,i);
    data = mean_ersp{i};  % alle Frequenzen bis 30Hz
    contourf(tf_times, tf_freqs, data, 50, 'LineStyle', 'none');
    clim([-2, 2])
    colormap('jet')
    xlabel('Time (ms)'); ylabel('Freq (Hz)');
    title(titles{i});
    colorbar;
end
sgtitle('Time-Frequency (<=30Hz) – SELF conditions');



% ======================
% OTHER Conditions (7–12) – TF Heatmaps bis 30Hz

figure;
for i = 7:12
    subplot(3,2,i-6);
     data = mean_ersp{i};  % alle Frequenzen bis 30Hz
    contourf(tf_times, tf_freqs, data, 50, 'LineStyle', 'none');
    clim([-2, 2])
    colormap('jet')
    xlabel('Time (ms)'); ylabel('Freq (Hz)');
    title(titles{i-6});
    colorbar;
end
sgtitle('Time-Frequency (<=30Hz) – OTHER conditions');

% ======================
% Theta (4–7Hz) Time Courses – SELF

figure; hold on;
colors = lines(6);

for i = 1:6
    theta_idx = tf_freqs >=4 & tf_freqs <=8;
    theta_mean = mean(mean_ersp{i}(theta_idx,:),1);  % mitteln über 4–7Hz
    plot(tf_times, theta_mean, 'LineWidth',1.5, 'Color',colors(i,:));
end
xlabel('Time (ms)');
ylabel('Theta Power (dB)');
legend(titles, 'Location','best');
title('Theta (4–7Hz) Time Courses – SELF conditions');
grid on;

% ======================
% Theta (4–7Hz) Time Courses – OTHER

figure; hold on;
colors = lines(6);

for i = 7:12
    theta_idx = tf_freqs >=4 & tf_freqs <=8;
    theta_mean = mean(mean_ersp{i}(theta_idx,:),1);  % mitteln über 4–7Hz
    plot(tf_times, theta_mean, 'LineWidth',1.5, 'Color',colors(i-6,:));
end
xlabel('Time (ms)');
ylabel('Theta Power (dB)');
legend(titles, 'Location','best');
title('Theta (4–7Hz) Time Courses – OTHER conditions');
grid on;


% Theta alle Lineplots in einer figure 

figure; hold on;
colors = lines(6);  % 6 Farben für die 6 Bedingungen

freq_band = [4 8];  % z.B. Theta
band_idx = tf_freqs >= freq_band(1) & tf_freqs <= freq_band(2);

yl = [-2 1]; % y-Achsenlimits nach deinem Datensatz
p1 = patch([50 100 100 50], [yl(1) yl(1) yl(2) yl(2)], [0.8 0.8 0.8], ...
    'FaceAlpha',0.3,'EdgeColor','none');
p2 = patch([275 375 375 275], [yl(1) yl(1) yl(2) yl(2)], [0.8 0.8 0.8], ...
    'FaceAlpha',0.3,'EdgeColor','none');

% Handles & Legenden
h = gobjects(1, n_cond); 
legend_labels = cell(1, n_cond);

for i = 1:n_cond
    band_mean = mean(mean_ersp{i}(band_idx,:), 1);

    if i <= 6
        linestyle = '--';  % Self: gestrichelt
        color = colors(i,:);
        legend_labels{i} = ['Self: ' titles{i}];
    else
        linestyle = '-';   % Other: durchgezogen
        color = colors(i-6,:);
        legend_labels{i} = ['Other: ' titles{i-6}];
    end

    h(i) = plot(tf_times, band_mean, 'LineWidth',1.5, ...
                'Color',color, 'LineStyle',linestyle);
end

xlabel('Time (ms)');
ylabel(sprintf('Power (%d–%d Hz, dB)', freq_band(1), freq_band(2)));
legend(h, legend_labels, 'Location', 'eastoutside', 'Interpreter','none');
title(sprintf('%d–%d Hz Power Time Courses – All 12 Conditions', freq_band(1), freq_band(2)));
ylim(yl);
grid on;

% ================
% Lineplots gezielt an elektroden FZ, FCz und Cz
channels = {'Fz', 'FCz', 'Cz'};
freq_band = [4 8];
band_idx = tf_freqs >= freq_band(1) & tf_freqs <= freq_band(2);
yl = [-2 1];

for c = 1:length(channels)
    chan = channels{c};
    chan_idx = find(strcmp({chanlocs.labels}, chan));  % Index des Kanals im ERSP-Datensatz
    
    figure('Name',['SELF/OTHER — ' chan]); hold on;
    colors = lines(6);
    h = gobjects(1, 12);
    legend_labels = cell(1,12);

    for i = 1:12
        % Zugriff direkt auf ERSP der gewünschten Bedingung und Kanal:
        ersp_data = squeeze(mean_ersp_raw{i});  % angenommen: ersp_raw{i} = [chan x freq x time]
        chan_data = squeeze(ersp_data(chan_idx,:,:));  % freq x time

        band_mean = mean(chan_data(band_idx,:), 1);  % mitteln über Frequenzband

        if i <= 6
            linestyle = '--';  % SELF
            color = colors(i,:);
            legend_labels{i} = ['SELF: ' titles{i}];
        else
            linestyle = '-';   % OTHER
            color = colors(i-6,:);
            legend_labels{i} = ['OTHER: ' titles{i-6}];
        end

        h(i) = plot(tf_times, band_mean, 'LineWidth', 1.5, ...
                    'Color', color, 'LineStyle', linestyle);
    end

    % Patches
    patch([50 100 100 50], [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], ...
          'FaceAlpha',0.3, 'EdgeColor','none');
    patch([275 375 375 275], [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], ...
          'FaceAlpha',0.3, 'EdgeColor','none');

    legend(h, legend_labels, 'Location','eastoutside', 'Interpreter','none');
    xlabel('Time (ms)');
    ylabel(sprintf('Power (%d–%d Hz, dB)', freq_band(1), freq_band(2)));
    title(sprintf('Theta (%d–%d Hz) – %s', freq_band(1), freq_band(2), chan));
    ylim(yl); grid on;
end

%%

% Alpha power 

% === Settings ===
frontal_roi = {'Fz', 'FCz', 'Cz'};
posterior_roi = {'Pz', 'POz', 'Oz'};
alpha_band = [8 12];

% === Helper function ===
compute_roi_mean = @(ersp, chan_roi) mean(ersp(chan_roi,:,:),1); % mitteln über Kanäle

% === Finde ROI Indizes ===
frontal_idx = find(ismember({chanlocs.labels}, frontal_roi));
posterior_idx = find(ismember({chanlocs.labels}, posterior_roi));

% === Frequenz & Zeit Indizes ===
alpha_idx = tf_freqs >= alpha_band(1) & tf_freqs <= alpha_band(2);

% === Farben für Plots ===
colors = lines(6);

% === SELF ===
figure;
hold on;

for i = 1:6
    data = mean_ersp{i}; % [freq x time]
    
    % Frontal
    frontal_alpha = squeeze(mean(data(alpha_idx,:),1)); % mitteln über Alpha
    plot(tf_times, frontal_alpha, 'LineWidth',1.5, 'Color',colors(i,:));
end
xlabel('Time (ms)');
ylabel('Alpha Power (dB)');
title('Alpha (8–12Hz) – SELF – Frontal Fz, FCz, Cz');
legend(titles(1:6), 'Location','best');
grid on;

% === OTHER ===
figure;
hold on;

for i = 1:6
    data = mean_ersp{i+6}; % [freq x time]
    
    % Frontal
    frontal_alpha = squeeze(mean(data(alpha_idx,:),1));
    plot(tf_times, frontal_alpha, 'LineWidth',1.5, 'Color',colors(i,:));
end
xlabel('Time (ms)');
ylabel('Alpha Power (dB)');
title('Alpha (8–12Hz) – OTHER – Frontal Fz, FCz, Cz');
legend(titles(1:6), 'Location','best');
grid on;

% === POSTERIOR ===
figure;
hold on;

for i = 1:6
    data = mean_ersp{i}; % SELF
    posterior_alpha = squeeze(mean(data(alpha_idx,:),1));
    plot(tf_times, posterior_alpha, 'LineWidth',1.5, 'Color',colors(i,:));
end
xlabel('Time (ms)');
ylabel('Alpha Power (dB)');
title('Alpha (8–12Hz) – SELF – Posterior Pz, POz, Oz');
legend(titles(1:6), 'Location','best');
grid on;

figure;
hold on;

for i = 1:6
    data = mean_ersp{i+6}; % OTHER
    posterior_alpha = squeeze(mean(data(alpha_idx,:),1));
    plot(tf_times, posterior_alpha, 'LineWidth',1.5, 'Color',colors(i,:));
end
xlabel('Time (ms)');
ylabel('Alpha Power (dB)');
title('Alpha (8–12Hz) – OTHER – Posterior Pz, POz, Oz');
legend(titles(1:6), 'Location','best');
grid on;

% Fertig
disp('Alle Plots erstellt.');

%% Grand averages 

% Get theta
theta=squeeze(mean(ersps(:,:,tf_freqs>=4&tf_freqs<=7,:), 3));

% ...at Fz
theta_fz=squeeze(theta(:,15,:));
theta_fcz=squeeze(theta(:,65,:));

% Average for incorrrect
theta_fz_incorrect = mean(theta_fz([1,2,7,8], :), 1);
theta_fz_flipped = mean(theta_fz([3,4,9,10], :), 1);

theta_fcz_incorrect = mean(theta_fcz([1,2,7,8], :), 1);
theta_fcz_flipped = mean(theta_fcz([3,4,9,10], :), 1);


figure()
plot(tf_times, theta_fz_incorrect, '-k', 'LineWidth', 2); hold on
plot(tf_times, theta_fz_flipped, ':k', 'LineWidth', 2)
plot(tf_times, theta_fcz_incorrect, '-r', 'LineWidth', 2)
plot(tf_times, theta_fcz_flipped, ':r', 'LineWidth', 2)

xlabel('Time (ms)')
ylabel('Theta Power (dB)')
title('Theta Power at Fz and FCz for Error Conditions')
legend({'Fz Incorrect', 'Fz Flipped', 'FCz Incorrect', 'FCz Flipped'}, ...
       'Location','best')
grid on


[~, peak_idx] = max(theta_fcz_incorrect(tf_times >= 0 & tf_times <= 100));
tmp = tf_times(tf_times >= 0 & tf_times <= 100); % ERN window 
peak_time_error = tmp(peak_idx);


[~, peak_idx] = max(theta_fcz_flipped(tf_times >= 300 & tf_times <= 600));
tmp = tf_times(tf_times >= 300 & tf_times <= 600); % FRN/P3 window
peak_time_flip = tmp(peak_idx);

figure()
plot(tf_times, theta_fz_incorrect, '-k', 'LineWidth', 2); hold on
plot(tf_times, theta_fz_flipped, ':k', 'LineWidth', 2)
plot(tf_times, theta_fcz_incorrect, '-r', 'LineWidth', 2)
plot(tf_times, theta_fcz_flipped, ':r', 'LineWidth', 2)

xlabel('Time (ms)')
ylabel('Theta Power (dB)')
title('Theta Power at Fz and FCz for Error Conditions')
legend({'Fz Incorrect', 'Fz Flipped', 'FCz Incorrect', 'FCz Flipped'}, ...
       'Location','best')
grid on

% Peak-Zeitfenster (±25 ms) grau hinterlegen
yL = ylim;  % aktuelle y-Achse
patch([peak_time_error-25 peak_time_error+25 peak_time_error+25 peak_time_error-25], ...
      [yL(1) yL(1) yL(2) yL(2)], ...
      [0.8 0.8 0.8], 'FaceAlpha',0.3, 'EdgeColor','none')

patch([peak_time_flip-25 peak_time_flip+25 peak_time_flip+25 peak_time_flip-25], ...
      [yL(1) yL(1) yL(2) yL(2)], ...
      [0.9 0.9 0.9], 'FaceAlpha',0.3, 'EdgeColor','none')

% vertikale Linien auf die Peaks selbst
xline(peak_time_error, 'y--', 'LineWidth', 1.5, 'DisplayName','Incorrect Peak');
xline(peak_time_flip,  'b--', 'LineWidth', 1.5, 'DisplayName','Flip Peak');

legend({'Fz Incorrect', 'Fz Flipped', 'FCz Incorrect', 'FCz Flipped', ...
        'Incorrect ROI', 'Flip ROI', 'Incorrect Peak', 'Flip Peak'}, 'Location','eastoutside')


% ROI-Frequenzbereich definieren (z.B. Theta oder Alpha)
freq_roi = [4 7];  % Beispiel: Theta
freq_idx = tf_freqs >= freq_roi(1) & tf_freqs <= freq_roi(2);

n_cond = 12;

% Array für die 12 Bedingungen (über Bedingungen und Frequenzen gemittelt)
all_cond = zeros(n_cond, length(tf_times));

for i = 1:n_cond
    tmp = mean_ersp{i}; % [freq x time]
    tmp_roi = mean(tmp(freq_idx, :), 1); % mitteln über Frequenzen
    all_cond(i, :) = tmp_roi; % speichern
end

% Grand Average über alle Bedingungen
grand_avg = mean(all_cond, 1); % [1 x time]


% Beispiel: GA aus deinem vorherigen Code
figure;
plot(tf_times, grand_avg, 'k', 'LineWidth', 2); hold on;
xlabel('Time (ms)');
ylabel('Power (dB)');
title('Grand Average Power (4–7 Hz)');
grid on;

% Search window 
search_window = [0 600];
search_idx = tf_times >= search_window(1) & tf_times <= search_window(2);

% Peak im Suchfenster finden
[peak_val, peak_idx_rel] = max(grand_avg(search_idx));


% Absoluter Index
search_times = tf_times(search_idx);
peak_time = search_times(peak_idx_rel);


% Define ±50 ms Zeitfenster um die gefundenen Peaks
roi_peak = [peak_time-50, peak_time+50];


% Patches zeichnen
yL = ylim;
patch([roi_peak fliplr(roi_peak)], [yL(1) yL(1) yL(2) yL(2)], ...
      [0.9 0.9 0.9], 'FaceAlpha',0.3,'EdgeColor','none');

% Textliche Ausgabe
fprintf('Peak ROI:   %.0f–%.0f ms (center at %.0f ms, max=%.2f dB)\n', ...
    roi_peak(1), roi_peak(2), peak_time, peak_val);

legend('Grand Average', 'Peak ROI');


%% Gezielte Avgerages für ANOVA

n_cond = 12;
cond_labels = { ...
   'Self_Falsch_Low', 'Self_Falsch_High', ...
   'Self_Flipped_Low', 'Self_Flipped_High', ...
   'Self_Korrekt_Low', 'Self_Korrekt_High', ...
   'Other_Falsch_Low', 'Other_Falsch_High', ...
   'Other_Flipped_Low', 'Other_Flipped_High', ...
   'Other_Korrekt_Low', 'Other_Korrekt_High' };

freq_roi = [4 7]; % Theta

% ROIs um Peaks (50ms Fenster)
error_roi = [peak_time_error-25, peak_time_error+25];
flip_roi  = [peak_time_flip-25, peak_time_flip+25];

chan_roi_error = {'Fz'};
chan_roi_flip  = {'FCz'};

tbl_error = table();
tbl_flip = table();

tbl_error.Subject = zeros(length(subject_list),1);
tbl_flip.Subject  = zeros(length(subject_list),1);

for s = 1:length(subject_list)
    subject = subject_list{s};
    fprintf('\n📄 Lade %s …\n', subject);

    load(fullfile(PATH_TF_DATA, [subject '_ersps_feedback_12cond.mat'])); % ersps, tf_freqs, tf_times, chanlocs

    if s == 1
        freq_idx = find(tf_freqs >= freq_roi(1) & tf_freqs <= freq_roi(2));
        time_idx_error = find(tf_times >= error_roi(1) & tf_times <= error_roi(2));
        time_idx_flip  = find(tf_times >= flip_roi(1)  & tf_times <= flip_roi(2));
        chan_idx_error = find(ismember({chanlocs.labels}, chan_roi_error));
        chan_idx_flip  = find(ismember({chanlocs.labels}, chan_roi_flip));
    end

    tbl_error.Subject(s) = str2double(subject(3:5));
    tbl_flip.Subject(s)  = str2double(subject(3:5));

    for c = 1:n_cond
        val_error = mean(ersps(c,chan_idx_error,freq_idx,time_idx_error), [2 3 4], 'omitnan');
        tbl_error.(cond_labels{c})(s,1) = val_error;

        val_flip = mean(ersps(c,chan_idx_flip,freq_idx,time_idx_flip), [2 3 4], 'omitnan');
        tbl_flip.(cond_labels{c})(s,1) = val_flip;
    end
end

outfile_error = fullfile(PATH_TF_DATA,'TF_ROI_JASP_ErrorPeak.csv');
outfile_flip  = fullfile(PATH_TF_DATA,'TF_ROI_JASP_FlipPeak.csv');

writetable(tbl_error, outfile_error);
writetable(tbl_flip,  outfile_flip);

fprintf('Error-Peak ROI-Tabelle gespeichert: %s\n', outfile_error);
fprintf('Flip-Peak ROI-Tabelle gespeichert: %s\n', outfile_flip);

disp(tbl_error);
disp(tbl_flip);


%% Topoplots
% Indizes für Theta & Zeitpunkte
freq_idx = find(tf_freqs >=4 & tf_freqs <=7);
[~, time_idx_error] = min(abs(tf_times - peak_time_error));
[~, time_idx_flip]  = min(abs(tf_times - peak_time_flip));

% Mittelung über Theta-Frequenz
mean_theta_at_error = squeeze(mean(ersps(:,:,freq_idx,time_idx_error),3));
mean_theta_at_flip  = squeeze(mean(ersps(:,:,freq_idx,time_idx_flip),3));

cond_labels = { ...
    'Self Falsch Low', 'Self Falsch High', ...
    'Self Flipped Low', 'Self Flipped High', ...
    'Self Korrekt Low', 'Self Korrekt High', ...
    'Other Falsch Low', 'Other Falsch High', ...
    'Other Flipped Low', 'Other Flipped High', ...
    'Other Korrekt Low', 'Other Korrekt High'};

%% === Error-Peak Topos ===
maxval_error = max(abs(mean_theta_at_error(:)));  % gleiche Skala für alle Subplots

figure('Name','Error-Peak Topoplots', ...
       'Position', [100 100 1400 800], ...
       'Color', 'w');
for cond = 1:12
    subplot(3,4,cond);
    data = mean_theta_at_error(cond,:);
    topoplot(data, chanlocs, ...
        'maplimits', [-maxval_error maxval_error], ...
        'electrodes', 'on');
    title(cond_labels{cond}, 'Interpreter', 'none', 'FontSize', 9);
end

sgtitle(sprintf('Theta Power Topos @ %.0f ms (Error-Peak)', peak_time_error), 'FontSize', 14);

cb = colorbar('Position', [0.92 0.1 0.02 0.8]);
cb.Label.String = 'Theta Power (dB)';
cb.FontSize = 10;

%% === Flip-Peak Topos ===
maxval_flip = max(abs(mean_theta_at_flip(:)));  

figure('Name','Flip-Peak Topoplots', ...
       'Position', [100 100 1400 800], ...
       'Color', 'w');
for cond = 1:12
    subplot(3,4,cond);
    data = mean_theta_at_flip(cond,:);
    topoplot(data, chanlocs, ...
        'maplimits', [-maxval_flip maxval_flip], ...
        'electrodes', 'on');
    title(cond_labels{cond}, 'Interpreter', 'none', 'FontSize', 9);
end

sgtitle(sprintf('Theta Power Topos @ %.0f ms (Flip-Peak)', peak_time_flip), 'FontSize', 14);

cb = colorbar('Position', [0.92 0.1 0.02 0.8]);
cb.Label.String = 'Theta Power (dB)';
cb.FontSize = 10;