%% Evaluate onset and offset tuning not sure 
% DEPENDENCIES:
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, February 2023


Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

% p_threshold = 0.05; % significance threshold for unit to be considered F0-sensitive

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

tuning = cell(20,4);
uCounter = 1;

load('all_peaks.mat')

totSpikeCount = [];

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp01/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    onset_tuning = zeros(length(units),length(stims),length(Flist))-1;
    offset_tuning = zeros(length(units),length(stims),length(Flist))-1;

    % for each unit
    for uu = 1:length(units)

        % if this unit wasn't active post-onset or offset, skip it
        if isempty(find(active(uu,:,:), 1))
            uCounter = uCounter + 1;
            continue
        end

        unitSpikes = Y(Y(:,3)==units(uu),:); % spikes for just this unit

        %%%%%%%%%%%%%%% ONSET TUNING %%%%%%%%%%%%%%%%%%
        if ~isempty(find(active(uu,:,1), 1))
            % for each stim type
            for ss = 1:length(stims)
                
                onset_response_spike_counts = zeros(length(repeats),length(Flist)); % initialize space to save the number of spikes evoked
                
                % for each F0
                for ff = 1:length(Flist)
    
                    stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); % unique name for combination of stim type and F0
        
                    if isempty(stimNum) % if this stim type and F0 combo wasn't presented
                        continue
                    end

                    if all_peaks(uCounter,ss,1) > 4
                        window = [0 0.1];
                    else
                        window = [0 0.05];
                    end
    
                    % for each presentation of this stim
                    for rr = 1:length(repeats)
    
                        onset_response_spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                        nSpikes = sum(onset_response_spikeIDXs);
                        totSpikeCount = [totSpikeCount; nSpikes];
                        onset_response_spike_counts(rr,ff) = nSpikes;
    
                    end % ends repeat loop
                end % ends F0 loop
            end % ends stim loop

            onset_tuning(uu,ss,:) = mean(onset_response_spike_counts);

        end % ends active onset if condition


        %%%%%%%%%%%%%%%% OFFSET TUNING %%%%%%%%%%%%%%%%%%%

        if ~isempty(find(active(uu,:,2), 1))
            % for each stim type
            for ss = 1:length(stims)
                
                offset_response_spike_counts = zeros(length(repeats),length(Flist));
                
                % for each F0
                for ff = 1:length(Flist)
    
                    stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); % unique name for combination of stim type and F0
        
                    if isempty(stimNum) % if this stim type and F0 combo wasn't presented
                        continue
                    end

                    if ap<10
                        if all_peaks(uCounter,ss,2) > 25
                            window = [0.3 0.4];
                        else
                            window = [0.3 0.35];
                        end
                    else
                        if all_peaks(uCounter,ss,2) > 18
                            window = [0.2 0.3];
                        else
                            window = [0.2 0.25];
                        end
                    end
    
                    % for each presentation of this stim
                    for rr = 1:length(repeats)
    
                        offset_response_spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                        nSpikes = sum(offset_response_spikeIDXs);
                        totSpikeCount = [totSpikeCount; nSpikes];
                        offset_response_spike_counts(rr,ff) = nSpikes;
    
                    end % ends repeat loop
                end % ends F0 loop

            end % ends stim loop

            offset_tuning(uu,ss,:) = mean(offset_response_spike_counts);

        end % ends active onset if condition       

        uCounter = uCounter + 1;

    end % ends unit loop

    tuning{ap,1} = Animals{ap};
    tuning{ap,2} = Pens{ap};
    tuning{ap,3} = onset_tuning;
    tuning{ap,4} = offset_tuning;

end % ends recording loop
