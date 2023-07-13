%% Plot click train tuning ordered by pitch salience to test given metric
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, April 2023

mat_struct = load('TNs_05_99.mat');
mat_cell = struct2cell(mat_struct);
units_by_rec = mat_cell{1};
clear mat_struct mat_cell

stims = {'CT0','CT5','CT10','CT20','CT40'};

nUnits = count_units(units_by_rec);
unit_counter = 1;
pitch_sensitivity = zeros(nUnits,1);
all_CT_tuning = cell(nUnits,1);

for pen = 1:length(units_by_rec)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/tmp02/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

    units = units_by_rec{pen,3};
    Flist = unique(F0);
    repeats = unique(Y(:,5));

    for uu = 1:size(units,1)

        unit = units(uu,1);

        unitSpikes = Y(Y(:,3)==unit,:); % get spikes for just this unit
        profile = zeros(length(stims),17);

        window = [0 .1];
        CT_tuning = cell(5,1);

         % go through each stim we want to plot
        for ss = 1:length(stims)
    
            nSpikes = zeros(length(repeats),length(Flist)); % allocate space to save spiking info for each trial
                
            for ff = 1:length(Flist) % for each frequency
        
                stimNum = find(strcmp(type,stims{ss}) & (F0==Flist(ff))); % unique name for combination of stim type and F0
    
                if isempty(stimNum) % if this stim type and fo combo wasn't presented
                    nSpikes(:,ff) = 0;
                    continue
                end
        
                % for each trial of this stim type
                for rr = 1:length(repeats)
                
                    spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                    nSpikes(rr,ff) = sum(spikeIDXs);
    
                end   
            end
    
            nSpikes = nSpikes ./ diff(window); % spikes per second
            
            CT_tuning{ss} = nSpikes;

        end % ends stim loop

        pitch_sensitivity(unit_counter) = estimate_pitch_salience_sensitivity_ANOVA(CT_tuning,0.01);
        all_CT_tuning{unit_counter} = CT_tuning;

        unit_counter = unit_counter + 1;
    end
end

%%

colors = colormap(hsv(length(stims)));

figure;
for uu = 1:nUnits

    if pitch_sensitivity(uu)
        figure(101);
    else
        figure(102);
    end

    nexttile; hold on;
    CT_tuning = all_CT_tuning{uu};

    for ss = 1:length(CT_tuning)
        responses = CT_tuning{ss};
        avg_responses = mean(responses);

        plot(1:17,avg_responses,'Color',colors(ss,:),'linewidth',1.5)
    end
    axis tight
    xticks([])
%     yticks([])
end