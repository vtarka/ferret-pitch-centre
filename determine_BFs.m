
%% Determine the best frequency for each condition for every unit
% BF defined as the F0 that evokes the largest mean spike rate
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, November 2022

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type); 
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));
    window = [0 0.15];

    % initialize some space to store the BFs
    BFs = zeros(length(units),length(stims));

    % for each unit
    for uu = 1:length(units)

        % if this unit was not frequency sensitive to any stim-type, skip it
        if isempty(find(sensitivity(uu,:), 1))
            sprintf('%s, %s, unit %d does not have any frequency selective neurons.',Animals{ap},Pens{ap},units(uu))
            continue
        end

        unitSpikes = Y(Y(:,3)==units(uu),:); % extract the spikes of just this unit

        % for each stim type
        for ss = 1:length(stims)
        
            % if the unit was frequency-sensitive to this condition, evaluate its best frequency
            if sensitivity(uu,ss)==1
                
                meanSpikes = zeros(length(Flist),1);

                % go through each F0 to find BF
                for ff = 1:length(Flist)

                    stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); % unique name for combination of stim type and F0
    
                    if isempty(stimNum) % if this stim type and fo combo wasn't presented
                        continue
                    end

                    nSpikes = zeros(length(repeats),1); % initialize space to count the number of spikes per trial repeat

                    % for each repeat
                    for rr = 1:length(repeats)

                        spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                        nSpikes(rr) = sum(spikeIDXs);

                    end % ends repeat loop

                    meanSpikes(ff) = mean(nSpikes); % average the spike rates across repeats

                end % ends looping through F0s
                
                [~,i] = max(meanSpikes); % find the index where the maximum spike rate occurs
                BFs(uu,ss) = i; % save the index of the best frequency

            end % ends if case
        end % ends loop through stims
    end % ends loop through units

    % UNCOMMENT BELOW TO SAVE THE BEST FREQUENCY VARIABLE IN THE SPIKING FILE
    % save(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat'],...
    % 'Y','type','F0','sensitivity','BFs')

end % ends loop through recordings