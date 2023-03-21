
%% Plot two tuning curves and their correlation starting from pairs with lowest correlation up to highest correlation
% Veronica Tarka
% veronica.tarka@dpag.ox.ac.uk
% January 2023

%%
Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

colors = colormap(jet(3));

max_rate_vs_rho = [];

for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));
    window = [0 0.15];

    stims_to_plot = {'CT0','tone'};

    % getting the indices for these stims
    stims_to_plot_IDX = zeros(length(stims_to_plot),1);
    for ss = 1:length(stims_to_plot)
        stims_to_plot_IDX(ss) = find(strcmp(stims,stims_to_plot{ss}));
    end

    unit_corrs = zeros(length(allUnits),2);

    for uu = 1:length(allUnits)
        unit = allUnits(uu);

        % check that the unit is actually F0-sensitive to one of these
        % stims
        if isempty(sensitivity(uu,stims_to_plot_IDX))
            continue
        end

        unitSpikes = Y(Y(:,3)==unit,:); % get the spikes of just this unit

        tuning = zeros(length(stims_to_plot),length(Flist));

        for ss = 1:length(stims_to_plot)
            
            nSpikes = zeros(length(repeats),length(Flist)); %creates an array of length repList by length Flist
            
            for ff = 1:length(Flist) % do the below for all the frequencies
        
                stimNum = find(strcmp(type,stims_to_plot{ss}) & (F0==Flist(ff))); 
    
                if isempty(stimNum) % if this stim type and fo combo wasn't presented
                    nSpikes(:,ff) = 0;
                    continue
                end
                
                % finds the stim label that corresponds to this stim type at
                % this particular F0
    
                % this stimNum will have been presented multiple times, 
                % so go through each presentation
        
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
        unit_corrs(uu,1) = unit;
        unit_corrs(uu,2) = rho;

    end % ends unit loop

    % get rid of all the units that weren't sensitive (flag=0)
    unit_corrs(:,unit_corrs(:,2)==0) = [];
    
    % order the correlations so we plot the lowest ones first
    ordered_unit_corrs = sortrows(unit_corrs,2);

    % figure out how many figures we'll need to show all these units
    nfigs = ceil(length(unit_corrs)/16);

    for fig = 1:nfigs
      
        figure('Position',[700 400 2100 1500]) 

        for sp = 1:16
            subplot(4,4,sp)

            if sp+16*(fig-1)>length(unit_corrs)
                continue
            end

            % get the tuning curves for these guys again
            unit = ordered_unit_corrs(sp+16*(fig-1),1);

            unitSpikes = Y(Y(:,3)==unit,:); % get the spikes of just this unit

            tuning = zeros(length(stims_to_plot),length(Flist));
    
            maxRate = 0;
            for ss = 1:length(stims_to_plot)
                
                nSpikes = zeros(length(repeats),length(Flist)); %creates an array of length repList by length Flist
                
                for ff = 1:length(Flist) % do the below for all the frequencies
            
                    stimNum = find(strcmp(type,stims_to_plot{ss}) & (F0==Flist(ff))); 
        
                    if isempty(stimNum) % if this stim type and fo combo wasn't presented
                        nSpikes(:,ff) = 0;
                        continue
                    end
                    
                    % finds the stim label that corresponds to this stim type at
                    % this particular F0
        
                    % this stimNum will have been presented multiple times, 
                    % so go through each presentation
            
                    for rr = 1:length(repeats)
                    
                        spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                        nSpikes(rr,ff) = sum(spikeIDXs);
        
                    end   
                end % ends F0 loop
        
                nSpikes = nSpikes ./ diff(window); % spikes per second

                meanSpikes = mean(nSpikes);
                if max(meanSpikes) > maxRate
                    maxRate = max(meanSpikes);
                end

                h = shadedErrorBar(1:17,meanSpikes,ste(nSpikes),{'Color',colors(ss,:)},1);
                h.mainLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
                h.edge(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
                h.edge(2).Annotation.LegendInformation.IconDisplayStyle = 'off';
                hold on
                if BFs(ss)~= 0
                    p = scatter(BFs(ss),meanSpikes(BFs(ss)),200,'MarkerFaceColor',colors(ss,:));
                    p.Annotation.LegendInformation.IconDisplayStyle = 'off';
                end
                xticks([])

                if sp==1
                    xticks(1:17)
                    xticklabels(num2str(Flist))
                    xlabel('F0')
                    ylabel('Evoked Firing Rate (spike/sec)'); % y axis label is the firing rate
                end
                title(sprintf('Rho = %.2f',ordered_unit_corrs(sp+16*(fig-1),2)))
                set(gca,'fontsize',16)
                axis tight
                hold on
    
            end % end stim loop

            max_rate_vs_rho = [max_rate_vs_rho; [ordered_unit_corrs(sp+16*(fig-1),2) maxRate]];
        end
        pause
    end
end % ends recording loop