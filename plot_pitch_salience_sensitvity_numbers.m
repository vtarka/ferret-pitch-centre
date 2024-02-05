
%% Plot CT0 tuning correlation with low and high tuning, pitch salience is the color map
% DEPENDENCIES:
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, May 2023

figure;

file_names = {'HN_units_new_05','TN_units_new_05','PN_units_new_05'};
% file_names = {'tone_HNs','tone_TNs','tone_PNs'};

stims = {'low','high','CT0','CT5','CT10','CT20','CT40'};


plot_titles = {'Harmonicity Neurons','Temporal Neurons','Pitch Neurons'};
% plot_titles = {'','',''};

window = [0 0.1];

max_sensitivity = 0;
min_sensitivity = 0;

colors = [0 0 1; 1 0 0; 0 0 0];

sensitivity_numbers = zeros(length(file_names),2);

for fi = 1:length(file_names)

    mat_struct = load(file_names{fi});
    mat_cell = struct2cell(mat_struct);
    units_by_rec = mat_cell{1};
    
    stims = {'CT0','CT5','CT10','CT20','CT40'};
    nUnits = count_units(units_by_rec);
    locs = [1 1 2 1 5 3 5 4 3 3 1 2 3 2 1 1 1 3 4 3];
    loc_list = zeros(nUnits,1);

    unit_counter = 1;
    pitch_sensitivity = zeros(nUnits,1);
    all_CT_tuning = cell(nUnits,1);
    
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
    
            pitch_sensitivity(unit_counter) = estimate_pitch_salience_sensitivity_ANOVA(CT_tuning,0.05);

            if pitch_sensitivity(unit_counter)
                loc_list(unit_counter) = locs(pen);
            end
            
            all_CT_tuning{unit_counter} = CT_tuning;
    
            unit_counter = unit_counter + 1;
        end
    end

    sensitivity_numbers(fi,1) = sum(pitch_sensitivity);
    sensitivity_numbers(fi,2) = count_units(units_by_rec) - sum(pitch_sensitivity);

    loc_list(loc_list==0) = [];
    figure;
    C = categorical(loc_list,1:5,{'low A1','high A1','low AAF','high AAF','PPF'});
    histogram(C,'FaceColor','k')
    ylabel('Number of Units')

    ylim([0 20])
    set(gca,'fontsize',18)

end % ends file loop

figure;

b = bar(sensitivity_numbers);
b(2).FaceColor = [0 0 0];
legend('Sensitive','Not Sensitive')
set(gca','fontsize',24)
xticklabels({'HNs','TNs','PNs'})
title('Pitch Salience Sensitivity')
ylabel('Neurons')


