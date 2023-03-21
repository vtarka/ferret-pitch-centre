
%% Veronica Tarka
% veronica.tarka@dpag.ox.ac.uk
% January 2023

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

totalPN_count = 0;
PN_CFs = [];
PN_units = cell(length(Animals),4);
PN_counter = 0;

for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));
    window = [0 0.15];

    stims_to_plot = {'CT0','tone'};

    % getting the indices for these stims
    stims_to_plot_IDX = zeros(length(stims_to_plot),1);
    for ss = 1:length(stims_to_plot)
        stims_to_plot_IDX(ss) = find(strcmp(stims,stims_to_plot{ss}));
    end

%     unit_corrs = zeros(length(allUnits),2);

    PN_unit_list = [];

    for uu = 1:length(units)
        unit = units(uu);

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
        if rho > 0.6
            PN_unit_list = [PN_unit_list; unit];
        end

    end % ends unit loop

    PN_units{ap,1} = Animals{ap};
    PN_units{ap,2} = Pens{ap};
    PN_units{ap,3} = PN_unit_list;

    PN_counter = PN_counter + length(PN_unit_list);
   
end % ends recording loop




%% Plot the tuning curves for the PNs we got from the above section

figure('Position',[1900 500 1800 1200])
sp = 1;
colors = colormap(jet(2));

for ap = 1:length(PN_units)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' PN_units{ap,1} '/tmp/Spikes_' PN_units{ap,1} '_' PN_units{ap,2} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));
    PNUnits = PN_units{ap,3};
    [~,PNUnit_IDXs] = ismember(PNUnits,allUnits);

    stims_to_plot = {'CT0','tone'};

    window = [0 0.15];

    for uu = 1:length(PNUnits)
        unit = PNUnits(uu);
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
    
            nSpikes = nSpikes ./ diff(window); % spikes per second

            subplot(4,4,sp)
            meanSpikes = mean(nSpikes);
            h = shadedErrorBar(1:17,meanSpikes,ste(nSpikes),{'Color',colors(ss,:)},1);
            h.mainLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
            h.edge(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
            h.edge(2).Annotation.LegendInformation.IconDisplayStyle = 'off';
            hold on
%             if BFs(ss)~= 0
%                 p = scatter(BFs(ss),meanSpikes(BFs(ss)),200,'MarkerFaceColor',colors(ss,:));
%                 p.Annotation.LegendInformation.IconDisplayStyle = 'off';
%             end
            xticks([])

            if sp==1
                xticks(1:17)
                xticklabels(num2str(Flist))
                xlabel('F0')
                ylabel('Evoked Firing Rate (spike/sec)'); % y axis label is the firing rate
            end
%             title(sprintf('Rho = %.2f',ordered_unit_corrs(sp+16*(fig-1),2)))
            set(gca,'fontsize',16)
            axis tight



        end % end stim loop

        if sp==16
            sp=1;
            figure('Position',[1900 500 1800 1200])
        else
            sp = sp+1;
        end
       
    end
end