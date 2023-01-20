
%% Look for pitch neurons as defined by a strong correlation between pure tone tuning and allHarm tuning
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, January 2023

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

rho_threshold = -1;
plot_yn = 'n'; % y = include plots of the unit's tuning, n = skip the plots

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

totalPN_count = 0;
test = 0;
PN_units = cell(length(Animals),3); % allocate space to save the units we find as harmonicity neurons
rhos = [];

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    stims_to_plot = {'allHarm','tone'};

    window = [0 0.15]; % response window

    PN_unit_list = []; % to keep track of HNs we find
    for uu = 1:length(units) % for each unit
        unit = units(uu);

        % if the unit is F0-sensitive to CT0 and low, but not to high, look at its correlation
        if sensitivity(uu,8)==1 && sensitivity(uu,13)==1

            test = test + 1;

            unitSpikes = Y(Y(:,3)==unit,:); % get the spikes of just this unit

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
                totalPN_count = totalPN_count + 1;
                PN_unit_list = [PN_unit_list, unit];
    
                if strcmp(plot_yn,'y')
    
                    stim_to_plot = {'allHarm','tone'};
                    BFs_to_plot = BFs(uu,[8 13]);
                    plot_tuning_by_cond(Y,type,F0,units(uu),stim_to_plot,BFs_to_plot,Animals{ap},Pens{ap});
    
                    pause
    
                end   
            end
        end
    end % ends the unit loop

    % save all the units we found for this penetration
    PN_units{ap,1} = Animals{ap};
    PN_units{ap,2} = Pens{ap};
    PN_units{ap,3} = PN_unit_list;

end % ends recording loop


%% Define it as a good correlation between low, high, and CT0

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

rho_threshold = -1;
plot_yn = 'n'; % y = include plots of the unit's tuning, n = skip the plots

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

totalPN_count = 0;
test = 0;
PN_units = cell(length(Animals),3); % allocate space to save the units we find as harmonicity neurons
rhos_low = [];
rhos_high = [];

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    stims_to_plot = {'low','high','CT0'};

    window = [0 0.15]; % response window

    PN_unit_list = []; % to keep track of HNs we find
    for uu = 1:length(units) % for each unit
        unit = units(uu);

        % if the unit is F0-sensitive to CT0 and low, but not to high, look at its correlation
        if sensitivity(uu,1)==1 && sensitivity(uu,10)==1 && sensitivity(uu,11)==1

            test = test + 1;

            unitSpikes = Y(Y(:,3)==unit,:); % get the spikes of just this unit

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
            rho_low = corr(tuning(1,:)',tuning(3,:)');
            rhos_low = [rhos_low; rho_low];

            rho_high = corr(tuning(2,:)',tuning(3,:)');
            rhos_high = [rhos_high; rho_high];

            if rho_low>rho_threshold && rho_high>rho_threshold
                totalPN_count = totalPN_count + 1;
                PN_unit_list = [PN_unit_list, unit];
    
                if strcmp(plot_yn,'y')
    
                    stim_to_plot = {'low','high','CT0'};
                    BFs_to_plot = BFs(uu,[11 10 1]);
                    plot_tuning_by_cond(Y,type,F0,units(uu),stim_to_plot,BFs_to_plot,Animals{ap},Pens{ap});
    
                    pause
    
                end   
            end
        end
    end % ends the unit loop

    % save all the units we found for this penetration
    PN_units{ap,1} = Animals{ap};
    PN_units{ap,2} = Pens{ap};
    PN_units{ap,3} = PN_unit_list;

end % ends recording loop