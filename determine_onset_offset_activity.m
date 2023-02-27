
%% Go through each stimulus type and determine whether the neuron is sensitive
% to variations in F0 using a one-way ANOVA
% DEPENDENCE: plot_tuning_by_cond
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, February 2023


Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

p_threshold = 0.05; % significance threshold for unit to be considered F0-sensitive

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    onset_window = [0 0.1];
    
    offset_window = [0.2 0.3];
    offset_window_N = [0.3 0.4];

    null_window = [0.6 0.7];
    null_window_N = [0.5 0.6];


    % want to build an nUnits by nConditions x 2 logical array (1 if active)
    active = zeros(length(units),length(stims),2);

    % for each unit
    for uu = 1:length(units)

        unitSpikes = Y(Y(:,3)==units(uu),:); % spikes for just this unit

        % for each stim type
        for ss = 1:length(stims)

            %%%%%%%%%%%%%%%% ONSET ACTIVITY %%%%%%%%%%%%%%%

            onset_response_spike_counts = zeros(length(repeats),length(Flist)); % initialize space to save the number of spikes evoked
            offset_response_spike_counts = zeros(length(repeats),length(Flist));
            null_spike_counts = zeros(length(repeats),length(Flist));

            % for each F0
            for ff = 1:length(Flist)

                stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); % unique name for combination of stim type and F0
    
                if isempty(stimNum) % if this stim type and F0 combo wasn't presented
                    continue
                end

                % for each presentation of this stim
                for rr = 1:length(repeats)
                
                    onset_response_spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr)...
                        & unitSpikes(:,2)>onset_window(1) & unitSpikes(:,2)<onset_window(2);
                    nSpikes = sum(onset_response_spikeIDXs);
                    onset_response_spike_counts(rr,ff) = nSpikes;

                    if strcmp(Animals{ap},'Noah')
                        null_spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr)...
                            & unitSpikes(:,2)>null_window_N(1) & unitSpikes(:,2)<null_window_N(2);
                        nSpikes = sum(null_spikeIDXs);
                        null_spike_counts(rr,ff) = nSpikes;

                        offset_response_spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr)...
                            & unitSpikes(:,2)>offset_window_N(1) & unitSpikes(:,2)<offset_window_N(2);
                        nSpikes = sum(offset_response_spikeIDXs);
                        offset_response_spike_counts(rr,ff) = nSpikes;
                    else
                        null_spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr)...
                            & unitSpikes(:,2)>null_window(1) & unitSpikes(:,2)<null_window(2);
                        nSpikes = sum(null_spikeIDXs);
                        null_spike_counts(rr,ff) = nSpikes;

                        offset_response_spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr)...
                            & unitSpikes(:,2)>offset_window(1) & unitSpikes(:,2)<offset_window(2);
                        nSpikes = sum(offset_response_spikeIDXs);
                        offset_response_spike_counts(rr,ff) = nSpikes;
                    end
                end 
            end

            onset_response_vector = reshape(onset_response_spike_counts,1,[]);
            offset_response_vector = reshape(offset_response_spike_counts,1,[]);
            null_vector = reshape(null_spike_counts,1,[]);

            % test whether this unit was significantly active in response to onset
            [hOn,~] = ttest2(onset_response_vector,null_vector,'alpha',p_threshold);
            [hOff,~] = ttest2(offset_response_vector,null_vector,'alpha',p_threshold);
            
            % if the p value is below our threshold, save the unit as F0-sensitive
            if hOn==1
                active(uu,ss,1) = 1;
            end 

            if hOff == 1
                active(uu,ss,2) = 1;
            end

        end % ends the stim loop  
    end % ends the unit loop

    % UNCOMMENT BELOW TO SAVE THE SENSITIVITY VARIABLE IN THE SPIKING FILE
    save(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp01/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat'],...
       'Y','type','F0','active')

end % ends the file-loading loop