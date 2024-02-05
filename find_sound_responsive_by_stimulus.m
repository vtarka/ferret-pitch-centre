% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, May 2023

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

sound_responsive = 0;
not_sound_responsive = 0;
alpha = 0.05;

responsive_units_per_pen = zeros(length(Animals),2);

for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp02/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type); 
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    if ap < 9
        pre_window = [0.47 0.57];
        post_window = [0.3 0.4];
    else
        pre_window = [0.65 0.75];
        post_window = [0.2 0.3];
    end

%     post_window = [0 0.1];

    % create variable to save whether each unit is sound responsive
    responsive_by_stim = zeros(length(units),length(stims));

    % for each unit
    for uu = 1:length(units)

        unitSpikes = Y(Y(:,3)==units(uu),:); % extract the spikes of just this unit

        % for each stim type
        for ss = 1:length(stims)

            pre_onset_spikes = [];
            post_onset_spikes = [];

            % go through each F0 to find BF
            for ff = 1:length(Flist)

                stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); % unique name for combination of stim type and F0

                if isempty(stimNum) % if this stim type and fo combo wasn't presented
                    continue
                end

                % for each repeat
                for rr = 1:length(repeats)

                    spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>pre_window(1) & unitSpikes(:,2)<pre_window(2);
                    pre_onset_spikes(end+1) = sum(spikeIDXs);

                    spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>post_window(1) & unitSpikes(:,2)<post_window(2);
                    post_onset_spikes(end+1) = sum(spikeIDXs);

                end % ends repeat loop

            end % ends looping through F0s

    
            if ttest2(pre_onset_spikes,post_onset_spikes,alpha) == 1
                responsive_by_stim(uu,ss) = 1; % mark this unit as responsive
            end

        end % ends loop through stims

    end % ends loop through units

%     responsive(responsive==0) = 1;
    % UNCOMMENT BELOW TO SAVE THE BEST FREQUENCY VARIABLE IN THE SPIKING FILE
    save(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp02/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat'],...
        'Y','type','F0','responsive_by_stim','responsive')

end % ends loop through recordings


%% plot responsiveness by stimulus


