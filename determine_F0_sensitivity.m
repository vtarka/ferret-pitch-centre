
%% Determine whether cells are F0-sensitive to each stimulus presented
% DEPENDENCIES: plot_tuning_by_cond.m (if plot_yn == 'y')
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, March 2023

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

p_threshold = 0.05; % significance threshold for unit to be considered F0-sensitive
plot_yn = 'n'; % y = include plots of every stimulus the unit is F0-sensitive to, n = skip the plots

% %stimList:         'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %  # (Noah #)       1       2          3         4        5             6            7             8 (10)       9 (11)    10 (12)   11 (13)  12 (14)   13 (15)

not_sensitive_u = []; % keep track of units found to not be sensitive at all
uCounter = 1;

window = [0 0.1]; % windows in seconds to evaluate F0 sensitivity

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp02/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    % want to build an nUnits by nConditions logical array (1 if F0-sensitive)
    sensitivity = zeros(length(units),length(stims));

    % for each unit
    for uu = 1:length(units)

        unitSpikes = Y(Y(:,3)==units(uu),:); % spikes for just this unit

        found_sensitive = 0; % flag to keep track of whether this unit was sensitive to any stimulus

        % for each stim type
        for ss = 1:length(stims)

            spike_counts = zeros(length(repeats),length(Flist)); % initialize space to save the number of spikes evoked
            group_labels = cell(1,length(Flist)); % initialize space to label each 'group' for the ANOVA (each F0 is separate group)

            % for each F0
            for ff = 1:length(Flist)

                stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); % unique name for combination of stim type and F0
    
                if isempty(stimNum) % if this stim type and F0 combo wasn't presented
                    group_labels{ff} = num2str(Flist(ff));
                    continue
                end

                % for each presentation of this stim
                for rr = 1:length(repeats)
                
                    spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<=window(2);
                    nSpikes = sum(spikeIDXs);

                    spike_counts(rr,ff) = nSpikes ./ diff(window); % convert to spike rate and save
    
                end 

                group_labels{ff} = num2str(Flist(ff)); % create group labels for MATLAB ANOVA function

            end

            % If this stimulus is low harmonics or masked low harmonics, the 1st two pitches were not presented
            if strcmp(stims{ss},'low') || strcmp(stims{ss},'F0MaskLow')
                spike_counts(:,1:2) = [];
                group_labels(1:2) = [];
            end

            % test whether F0 significantly modulates spike rate
            [p,tbl,stats] = anova1(spike_counts,group_labels,'off');
            
            % if the p value is below our threshold, save the unit as F0-sensitive
            if p < p_threshold
                sensitivity(uu,ss) = 1;
                found_sensitive = 1;
            end 

        end % ends the stim loop  

        if found_sensitive == 0 
            not_sensitive_u(end+1) = uCounter;
        end

        if strcmp(plot_yn,'y')
            if ~isempty(find(sensitivity(uu,:), 1))

                stim_to_plot = {};
                for ss = 1:length(stims)
                
                    if sensitivity(uu,ss)==1
                        stim_to_plot{end+1} = stims{ss};
                    end

                end

                clf
                plot_tuning_by_cond(Y,type,F0,units(uu),stim_to_plot,zeros(length(stim_to_plot),1),Animals{ap},Pens{ap},window);

                pause

            end % ends no sensitivity if case
        end % ends plot if case

        uCounter = uCounter + 1;

    end % ends the unit loop

    % UNCOMMENT BELOW TO SAVE THE SENSITIVITY VARIABLE IN THE SPIKING FILE
    save(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/p05/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat'],...
       'Y','type','F0','sensitivity')

end % ends the file-loading loop