%% Find harmonicity neurons by F0-sensitivity to high, low, and CT0
% DEPENDENCIES: plot_tuning_by_cond.m (if plot_yn == 'y')
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, March 2023

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};


% %stimList:         'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %  # (Noah #)       1       2          3         4        5             6            7             8 (10)       9 (11)    10 (12)   11 (13)  12 (14)   13 (15)


stims_to_plot = {'CT0','alt','high','low','tone'};

onset_windows = [0:5:95; 5:5:100]'/1000; % in ms

colors = colormap(hsv(length(stims_to_plot)));

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp02/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    if ap<9
        offset_windows = [300:5:395; 305:5:400]'/1000;
    else
        offset_window = [200:5:295; 205:5:300]'/1000;
    end

    stims = unique(type);
    units = unique(Y(:,3));
    Flist = unique(F0);
    repeats = unique(Y(:,5));

    for uu = 1:length(units) % for each unit
        unit = units(uu);

        if responsive(uu)==0
            continue
        end

        unitSpikes = Y(Y(:,3)==units(uu),:);
        onset_PSTHs = zeros(length(stims_to_plot),length(onset_windows));
        offset_PSTHs = zeros(length(stims_to_plot),length(onset_windows));

        for ss = 1:length(stims_to_plot)

            onset_spike_counts = zeros(length(Flist)*length(repeats),length(onset_windows));
            offset_spike_counts = zeros(length(Flist)*length(repeats),length(offset_windows));
            trial_counter = 1;
                        
            for ff = 1:length(Flist)

                stimNum = find(strcmp(type,stims_to_plot(ss)) & (F0==Flist(ff))); % unique name for combination of stim type and F0
    
                if isempty(stimNum) % if this stim type and F0 combo wasn't presented
                    continue
                end

                
                % for each presentation of this stim
                for rr = 1:length(repeats)

                    for ww = 1:length(onset_windows)
                
                        onset_spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>onset_windows(ww,1) & unitSpikes(:,2)<=onset_windows(ww,2);
                        offset_spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>offset_windows(ww,1) & unitSpikes(:,2)<=offset_windows(ww,2);
                        onset_nSpikes = sum(onset_spikeIDXs);
                        offset_nSpikes = sum(offset_spikeIDXs);
    
                        onset_spike_counts(trial_counter,ww) = onset_nSpikes/0.005; % convert to spike rate and save
                        offset_spike_counts(trial_counter,ww) = offset_nSpikes/0.005;
                    end

                    trial_counter = trial_counter + 1;   

                end % ends repeat loop
            end % ends frequency loop

            onset_PSTHs(ss,:) = mean(onset_spike_counts);
            offset_PSTHs(ss,:) = mean(offset_spike_counts);

            subplot(2,1,1); hold on
            plot(1:length(onset_windows),mean(onset_spike_counts),'color',colors(ss,:),'linewidth',2)

            xticks(0:5:20)
            xticklabels(0:5*5:20*5)

            subplot(2,1,2); hold on
            plot(1:length(offset_windows),mean(offset_spike_counts),'color',colors(ss,:),'linewidth',2)
            xticks(0:5:20)
            xticklabels(0:5*5:20*5)

        end % ends stim loop

        max_onset_rate = max(max(onset_PSTHs));
        max_offset_rate = max(max(offset_PSTHs));
        ymax = max(max_onset_rate,max_offset_rate);

        subplot(2,1,2)
        legend(stims_to_plot)
        ylim([0 ymax+(ymax*0.1)])

        subplot(2,1,1)
        ylim([0 ymax+(ymax*0.1)])


        pause;
        clf

    end % ends unit loop

end % ends recording loop
