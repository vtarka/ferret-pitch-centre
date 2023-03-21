
%% Determine the best frequency for each condition for every unit
% BF defined as the F0 that evokes the largest mean spike rate


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

    BFs = zeros(length(units),length(stims));

    % for each unit
    for uu = 1:length(units)

        if isempty(find(sensitivity(uu,:), 1))
            sprintf('%s, %s, unit %d does not have any frequency selective neurons.',Animals{ap},Pens{ap},units(uu))
            continue
        end

        unitSpikes = Y(Y(:,3)==units(uu),:);

        for ss = 1:length(stims)
        
            if sensitivity(uu,ss)==1
                
                meanSpikes = zeros(length(Flist),1);

                % go through each F0 to find BF
                for ff = 1:length(Flist)

                    stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); 
    
                    if isempty(stimNum) % if this stim type and fo combo wasn't presented
                        continue
                    end

                    nSpikes = zeros(length(repeats),1);

                    for rr = 1:length(repeats)
                        spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                        nSpikes(rr) = sum(spikeIDXs);

                    end % ends looping through stim repeats

                    meanSpikes(ff) = mean(nSpikes);

                end % ends looping through F0s
                
                [~,i] = max(meanSpikes);
                BFs(uu,ss) = i;

            end % ends if case
        end % ends loop through stims
    end % ends loop through units

    save(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat'],...
    'Y','type','F0','sensitivity','BFs')

end % ends loop through recordings

