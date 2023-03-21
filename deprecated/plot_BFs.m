

%% Veronica Tarka
% December 2022
% veronica.tarka@dpag.ox.ac.uk

%% Look at units that are frequency sensitive to both pure tones and CT0
% Plot tuning curves with a big marker over the best frequency

% List of all the recordings separating out the animal name and penetration
% intended to be looped through together (i.e. Animals{1},Pens{1} is a
% pair)
Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

% for every recording
for ap = 1:length(Animals)

%     if ap<12
%         continue
%     end

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);

    % stims to be included in the plot (choose from the stim list)
    stims_to_plot = {'low'};

    % getting the indices for these stims
    stims_to_plot_IDX = zeros(length(stims_to_plot),1);
    for ss = 1:length(stims_to_plot)
        stims_to_plot_IDX(ss) = find(strcmp(stims,stims_to_plot{ss}));
    end

    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));
    window = [0 0.15]; % window of time where a spike is considered a response

    % for each unit
    for uu = 1:length(units)

        % if this unit was not frequency-sensitive in response to any
        % condition, skip this unit
        if isempty(find(sensitivity(uu,stims_to_plot_IDX), 1))
            continue
        end

        BFs_to_plot = zeros(length(stims_to_plot),1);
        for ss = 1:length(stims_to_plot)
            BFs_to_plot(ss) = BFs(uu,stims_to_plot_IDX(ss));
        end

        plot_tuning_by_cond(Y,type,F0,units(uu),stims_to_plot,BFs_to_plot,Animals{ap},Pens{ap});
        pause

    end % ends unit loop
end % ends recording loop


%% Plot the stims eliciting frequency sensitivity topologically
% each stim is a node, edge between the two nodes if a unit is tuned at
% both

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

figure('Position',[1500 500 1800 1200])



for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    
    stims = unique(type);
%     stim_to_plot_IDX = [find(strcmp(stims,'low')), find(strcmp(stims,'high')), find(strcmp(stims,'CT0')),find(strcmp(stims,'tone'))];
%     stim_to_plot = {'low','high','CT0','tone'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));
    window = [0 0.15];


    % for each unit
    for uu = 1:length(units)

        if isempty(find(sensitivity(uu,:), 1))
            continue
        end

        a = zeros(length(stims));
        clf; sgtitle(sprintf('%s, %s unit %d',Animals{ap},Pens{ap},units(uu)))

        sensitive_stims = find(sensitivity(uu,:));

        if length(sensitive_stims)==1
            a(sensitive_stims,sensitive_stims) = a(sensitive_stims,sensitive_stims) + 1;
        else
            for ss = 1:length(sensitive_stims)
                for s = 1:length(sensitive_stims)
                    if sensitive_stims(ss) ~= sensitive_stims(s)
                        a(sensitive_stims(ss),sensitive_stims(s)) = a(sensitive_stims(ss),sensitive_stims(s)) + 1;
                    end
                end
            end
        end

    
        G = graph(a,stims);
        LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
        plot(G,'LineWidth',LWidths,'layout','circle')
        set(gca,'fontsize',20)
        pause

    end % ends unit loop

end % ends recording loop


%% Plot response profiles for single units
% Each plot will be an image, with the x-axis being F0 and each row being
% one condition type. The best frequency (if existing) will be highlighted
% in a red square. The color will indicate z-score.



Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

figure('Position',[2100 500 1600 1000])



for ap = 1:length(Animals)

    if ap<9
        continue
    end

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    
    stims = unique(type);
    stims_to_plot = {'CT0','CT5','CT10','CT20','CT40'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));
    window = [0 0.15];

    stim_to_plot_IDX = [];
    for ss = 1:length(stims_to_plot)
        stim_to_plot_IDX = [stim_to_plot_IDX; find(strcmp(stims,stims_to_plot{ss}))];
    end

    % for each unit
    for uu = 1:length(units)
        unit = units(uu); % get this unit's number

        % if this unit wasn't F0-sensitive for any condition, skip it
        if isempty(find(sensitivity(uu,:), 1))
            continue
        end

        % get our figure set up
        clf;
        sgtitle(sprintf('%s, %s unit # %d',Animals{ap},Pens{ap},unit))

        unitSpikes = Y(Y(:,3)==unit,:); % get the spikes of just this unit

        % make an empty array to store the z-scored response to each stim
        respProf = zeros(length(stims_to_plot),length(Flist));

        % get the best frequencies for all the stims we're plotting
        % if the unit wasn't sensitive to a particular stim, the BF will be
        % 0
        BFs_to_plot = [];
        for bf_IDX = 1:length(stim_to_plot_IDX)
            BFs_to_plot = [BFs_to_plot,BFs(uu,stim_to_plot_IDX(bf_IDX))];
        end % ends loop through best frequencies

        if isempty(find(BFs_to_plot))
            continue
        end

        for ss = 1:length(stims_to_plot)

            nSpikes = zeros(length(repeats),length(Flist));

            for ff = 1:length(Flist)
    
                stimNum = find(strcmp(type,stims_to_plot(ss)) & (F0==Flist(ff))); 
    
                if isempty(stimNum) % if this stim type and fo combo wasn't presented
                    continue
                end
            
                for rr = 1:length(repeats)

                    spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                    nSpikes(rr,ff) = sum(spikeIDXs);

                end % ends repeat loop
            end % ends F0 loop

            meanSpikes = mean(nSpikes);
            respProf(ss,:) = zscore(meanSpikes);

        end % ends stim loop

        colormap default;
        imagesc(respProf);  
        cb = colorbar;
        cb.Ticks = floor(min(min(respProf))):ceil(max(max(respProf)));
        set(gca,'YDir','reverse')
        xticks(1:17); xticklabels(Flist)
        yticks(1:length(stims_to_plot)); yticklabels(stims_to_plot)
        hold on; 
        colormap default;
        
        for bf = 1:length(BFs_to_plot)
            plot(BFs_to_plot(bf),bf,'r+', 'LineWidth',10,'MarkerSize',30)
        end % ends loop through BFs 

        plot_tuning_by_cond(Y,type,F0,unit,stims_to_plot,BFs(uu,stim_to_plot_IDX),Animals{ap},Pens{ap});
        figure(1); colormap default;
        pause
    

    end % ends unit loop
end % ends recording loop
