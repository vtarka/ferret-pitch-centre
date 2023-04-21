%% Plot how sensitive each temporal neuron is to pitch salience vs sensitivity to phase shifts (alt/rand)
% DEPENDENCIES: estimate_pitch_salience_sensitivity.m, estimate_phase_shift_sensitivity.m, count_units.m
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, April 2023

load('TNs_fNoah')
units_by_rec = TN_units;

stims = {'CT0','CT5','CT10','CT20','CT40','high','alt','rand'};

nUnits = count_units(units_by_rec);
phase_and_pitch_sensitivity = zeros(nUnits,2)-2;
unit_counter = 1;
real_profiles = zeros(nUnits,length(stims),17);

locs = [1 1 2 1 5 3 5 4 3 3 1 2 3 2 1 1 1 3 4 3];
loc_list = zeros(nUnits,1);

for pen = 1:length(units_by_rec)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/tmp02/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

    units = units_by_rec{pen,3};

    if ~isempty(units)
        units(units(:,2)==3,:) = [];
        units = unique(units(:,1));
    else
        continue
    end

    Flist = unique(F0);
    repeats = unique(Y(:,5));

    for uu = 1:size(units,1)

        unit = units(uu,1);

        unitSpikes = Y(Y(:,3)==unit,:); % get spikes for just this unit
        profile = zeros(length(stims),17);
% 
%         if units(uu,2)==1
%             window = [0 0.06];
%         elseif units(uu,2)==2
%             window = [0.06 0.15];
%         else
%             if pen<9
%                 window = [.3 .4];
%             else
%                 window = [.2 .3];
%             end
%         end

        window = [0 0.06];

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
            meanSpikes = mean(nSpikes); % average across repeats
            profile(ss,:) = meanSpikes;

            if strcmp('low',stims{ss})
                profile(ss,1:2) = nan;
            end

        end % ends stim loop

        phase_and_pitch_sensitivity(unit_counter,1) = estimate_pitch_sensitivity(profile(1:5,:));
        phase_and_pitch_sensitivity(unit_counter,2) = estimate_phase_shift_sensitivity(profile(6:8,:));
        real_profiles(unit_counter,:,:) = profile;

        loc_list(unit_counter) = locs(pen);
        unit_counter = unit_counter + 1;
    end
end

% not_phase_sensitive = find(phase_and_pitch_sensitivity(:,2)==-1);
% phase_sensitive = find(phase_and_pitch_sensitivity(:,2)==1);
% 
% figure;
% scatter(zeros(length(not_phase_sensitive),1),phase_and_pitch_sensitivity(not_phase_sensitive,1))
% 
% hold on;
% scatter(ones(length(phase_sensitive),1),phase_and_pitch_sensitivity(phase_sensitive,1));
% 
% figure;
% for pp = 1:length(phase_sensitive)
%     nexttile;
%     imagesc(squeeze(real_profiles(phase_sensitive(pp),1:5,:)))
% end