%% Plot the PSTHs for very responsive units to get a better sense of their latency and response duration
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, March 2023


Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

bin_size_ms = 5;
windows = [0:bin_size_ms:300; bin_size_ms:bin_size_ms:305]'/1000; % 10 ms windows to sum spiketimes over to build the PSTH

figure;

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp02/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    for uu = 1:length(units)

        unitSpikes = Y(Y(:,3)==units(uu),:);

        [r,c] = find(unitSpikes(:,2));
        all_trials_with_spikes = unique(unitSpikes(r,6));

        if length(find(unitSpikes(:,2)<0.3 & unitSpikes(:,2)>0)) > 250

            trials = unique(unitSpikes(:,6));
            nTrials = length(trials);
    
            nSpikes = zeros(nTrials,length(windows));
    
            for ww = 1:length(windows)
    
                win_start = windows(ww,1);
                win_end = windows(ww,2);
    
                [r,c] = find(unitSpikes(:,2) > win_start & unitSpikes(:,2) <= win_end);
                trials_with_spikes = unique(unitSpikes(r,6));
    
                for tt = 1:nTrials
    
                    if ismember(trials(tt),trials_with_spikes)
    
                        spikes = unitSpikes(:,2) > win_start & unitSpikes(:,2) <= win_end & unitSpikes(:,6) == trials(tt);
                        nSpikes(tt,ww) = sum(spikes);
                    end
    
                end % ends trial loop
            end % ends window loop

            avgPSTH = mean(nSpikes);

            clf
            plot(avgPSTH,'linewidth',3)
            axis tight
            xticks(0:5:60)
            xticklabels(0:5*5:60*5)
            pause


        end % ends if case
    end % ends unit loop
end % ends recording loop