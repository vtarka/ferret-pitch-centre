%% Find periodicity (temporal) neurons by responses within response windows defined by the cluster assignments
% Cluster 1: [0 80] 
% Cluster 2: [0 150]
% Cluster 4: [200 300] (unless Noah, then [300 400])
% DEPENDENCIES: 
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, March 2023

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

locs = [1 1 2 1 5 3 5 4 3 3 1 2 3 2 1 1 1 3 4 3]; 

Qualia = 'Good';

rho_threshold = -1;
plot_yn = 'n'; % y = include plots of the unit's tuning, n = skip the plots

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

totalTN_count = 0;
TN_units = cell(length(Animals),3); % allocate space to save the units we find as harmonicity neurons
rhos = [];
c_counts = [];
loc_counts = [];

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp01/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    stims_to_plot = {'CT0','high'};

    TN_unit_list = []; % to keep track of HNs we find
    for uu = 1:length(units) % for each unit
        unit = units(uu);

        cs = cluster_sensitivity;

        % if the unit is F0-sensitive to CT0 and low, but not to high, look at its correlation
%         if cluster_sensitivity(uu,1)==1 && cluster_sensitivity(uu,10)==1
        if cs(uu,1) && ~cs(uu,3) && ~cs(uu,4) && cs(uu,10)

            c_counts(end+1) = clusters(uu);
            loc_counts(end+1) = locs(ap);

            unitSpikes = Y(Y(:,3)==unit,:); % get the spikes of just this unit

            if clusters(uu) == 1
                window = [0 0.08];
            elseif clusters(uu) == 2 
                window = [0 0.15];
            elseif clusters(uu) == 4
                if ap < 9 
                    window = [0.3 0.4];
                else
                    window = [0.2 0.3];
                end
            end

            tuning = zeros(length(stims_to_plot),length(Flist));
    
            % for each stim
            for ss = 1:length(stims_to_plot)
                
                nSpikes = zeros(length(repeats),length(Flist)); % create an emtpy matrix to hold spiking in each trial
                
                % for each F0
                for ff = 1:length(Flist)
            
                    stimNum = find(strcmp(type,stims_to_plot{ss}) & (F0==Flist(ff)));  % unique name for combination of stim type and F0
        
                    if isempty(stimNum) % if this stim type and fo combo wasn't presented
                        nSpikes(:,ff) = 0;
                        continue
                    end

                    % for each presentation of this stim
                    for rr = 1:length(repeats)
                    
                        spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                        nSpikes(rr,ff) = sum(spikeIDXs);
        
                    end   
                end % ends F0 loop
        
                nSpikes = mean(nSpikes ./ diff(window)); % spikes per second
                tuning(ss,:) = nSpikes;
    
            end % end stim loop
    
            % now evaluate the correlation
            rho = corr(tuning(1,:)',tuning(2,:)');
            rhos = [rhos; rho];

            if rho>rho_threshold
                totalTN_count = totalTN_count + 1;
                TN_unit_list = [TN_unit_list, unit];
    
                if strcmp(plot_yn,'y')
    
                    clf
                    stim_to_plot = {'high','low','CT0'};
                    BFs_to_plot = [0 0 0 ]; %BFs(uu,[10 11 1]);
                    plot_tuning_by_cond(Y,type,F0,units(uu),stim_to_plot,BFs_to_plot,Animals{ap},Pens{ap});
                    sgtitle(sprintf('Cluster %d',clusters(uu)))
    
                    pause
    
                end   
            end
        end
    end % ends the unit loop

    % save all the units we found for this penetration
    TN_units{ap,1} = Animals{ap};
    TN_units{ap,2} = Pens{ap};
    TN_units{ap,3} = TN_unit_list;

end % ends recording loop


%% Check for tuning adjustments to high vs alt vs rand

for pen = 1:length(TN_units)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' TN_units{pen,1} '/tmp/Spikes_' TN_units{pen,1} '_' TN_units{pen,2} '_Good_Pitch.mat']);

    units = TN_units{pen,3};
    stim_to_plot = {'high','alt','rand'};
    BFs_to_plot = [0 0 0];

    for uu = 1:length(units)
    
        clf
        plot_tuning_by_cond(Y,type,F0,units(uu),stim_to_plot,BFs_to_plot,Animals{pen},Pens{pen});
        pause

    end

end


%% Check for tuning at click train jitters beyond 10%

for pen = 1:length(TN_units)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' TN_units{pen,1} '/tmp/Spikes_' TN_units{pen,1} '_' TN_units{pen,2} '_Good_Pitch.mat']);

    units = TN_units{pen,3};
    stim_to_plot = {'CT0','CT5','CT10','CT20','CT40'};
    BFs_to_plot = [0 0 0 0 0];

    for uu = 1:length(units)
    
        clf
        plot_tuning_by_cond(Y,type,F0,units(uu),stim_to_plot,BFs_to_plot,Animals{pen},Pens{pen});
        pause

    end

end
