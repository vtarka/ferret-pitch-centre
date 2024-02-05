%% Plot click train tuning ordered by pitch salience to test given metric
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, April 2023

mat_struct = load('TN_units_05_99_2.mat');
mat_cell = struct2cell(mat_struct);
units_by_rec = mat_cell{1};
clear mat_struct mat_cell

stims = {'high','alt','rand'};

nUnits = count_units(units_by_rec);
unit_counter = 1;
phase_sensitivity = zeros(nUnits,1);
all_phase_tuning = cell(nUnits,1);

for pen = 1:length(units_by_rec)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/final/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

    units = units_by_rec{pen,3};
    Flist = unique(F0);
    repeats = unique(Y(:,5));

    for uu = 1:size(units,1)

        unit = units(uu,1);

        unitSpikes = Y(Y(:,3)==unit,:); % get spikes for just this unit
        profile = zeros(length(stims),17);

        window = [0 .1];
        phase_tuning = cell(length(stims),1);

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
            
            phase_tuning{ss} = nSpikes;

            if isempty(phase_tuning{ss})
                a = 1;
            end

        end % ends stim loop

        phase_sensitivity(unit_counter) = estimate_phase_shift_sensitivity_ANOVA(phase_tuning,0.01);
        all_phase_tuning{unit_counter} = phase_tuning;

        unit_counter = unit_counter + 1;
    end
end

%%

colors = colormap(hsv(length(stims)));

figure;
for uu = 1:nUnits

    if phase_sensitivity(uu)
        figure(101);
    else
        figure(102);
    end

    nexttile; hold on;
    phase_tuning = all_phase_tuning{uu};

    for ss = 1:length(phase_tuning)
        responses = phase_tuning{ss};
        avg_responses = mean(responses);

        plot(1:17,avg_responses,'Color',colors(ss,:),'linewidth',1.5)
    end
    axis tight
    xticks([])
%     yticks([])
end