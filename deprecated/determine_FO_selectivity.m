
%% Go through each condition type and determine whether the neuron is sensitive
% to variations in F0 using a one-way ANOVA


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

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));
    window = [0 0.15];

    % want to build an nUnits by nConditions logical array (1 if
    % F0-sensitive)
    sensitivity = zeros(length(units),length(stims));

    for uu = 1:length(units)

        unitSpikes = Y(Y(:,3)==units(uu),:);

        % for each stim type
        for ss = 1:length(stims)

            spike_counts = zeros(length(repeats),length(Flist));
            group_labels = cell(1,length(Flist));

            % for each F0
            for ff = 1:length(Flist)

                stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); 
    
                if isempty(stimNum) % if this stim type and F0 combo wasn't presented
                    group_labels{ff} = num2str(Flist(ff));
                    continue
                end

                for rr = 1:length(repeats)
                
                    spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                    nSpikes = sum(spikeIDXs);

                    spike_counts(rr,ff) = nSpikes ./ diff(window);
    
                end 

                group_labels{ff} = num2str(Flist(ff));

            end

            %%% DO THE ANOVA HERE
            [p,tbl,stats] = anova1(spike_counts,group_labels,'off');
            
            if p<0.05
                sensitivity(uu,ss) = 1;
            end 

        end % ends the stim loop      
    end % ends the unit loop

    save(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat'],...
        'Y','type','F0','sensitivity')

end % ends the file-loading loop





%% Investigate individual units

% load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

stims = unique(type);
Flist = unique(F0);
repeats = unique(Y(:,5));
units = unique(Y(:,3));
window = [0 0.15];

for uu = 1:length(units)

    if ~isempty(find(sensitivity(uu,:), 1))
        sprintf('\nUnit: %d:\n',units(uu))

        stim_to_plot = {};
        for ss = 1:length(stims)
        
            if sensitivity(uu,ss)==1
                disp(stims{ss})
                stim_to_plot{end+1} = stims{ss};

            end

        end

        plot_tuning_by_cond(Y,type,F0,units(uu),stim_to_plot);
        pause

    end
end


















