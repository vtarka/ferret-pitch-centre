%% Plot click train tuning ordered by pitch salience to test given metric
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, April 2023



load('TNs_fNoah')
units_by_rec = TN_units;

stims = {'CT0','CT5','CT10','CT20','CT40'};

nUnits = count_units(units_by_rec);
unit_counter = 1;
CT_profiles = zeros(nUnits,length(stims),17);
pitch_sensitivity = zeros(nUnits,1);

for pen = 1:length(units_by_rec)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/tmp02/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

    units = units_by_rec{pen,3};
    Flist = unique(F0);
    repeats = unique(Y(:,5));

    for uu = 1:size(units,1)

        unit = units(uu,1);

        unitSpikes = Y(Y(:,3)==unit,:); % get spikes for just this unit
        profile = zeros(length(stims),17);

        if units(uu,2)==1
            window = [0 0.06];
        elseif units(uu,2)==2
            window = [0.06 0.15];
        else
            if pen<9
                window = [.3 .4];
            else
                window = [.2 .3];
            end
        end

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

        pitch_sensitivity(unit_counter) = estimate_pitch_sensitivity(profile);
        CT_profiles(unit_counter,:,:) = profile;

        unit_counter = unit_counter + 1;
    end
end

[B,I] = sort(pitch_sensitivity);
colors = colormap(hsv(length(stims)));

figure;
for uu = 1:nUnits
    nexttile; hold on;
    for ss = 1:length(stims)
        plot(1:17,squeeze(CT_profiles(I(uu),ss,:)),'Color',colors(ss,:),'linewidth',1.5)
    end
    axis tight
    xticks([])
    yticks([])
    title(sprintf('%.2f',B(uu)))
end



function sensitivity = estimate_pitch_sensitivity(CT_tuning)

% CT_tuning should be 5 x 17 where each row is the tuning curve for one of
% the CT stimuli IN ORDER so that CT_tuning(1,:) = tuning for CT0,
% CT_tuning(2,:) = tuning for CT5, etc up to CT40

% max_spike_rate = max(max(CT_tuning));
% CT_tuning = CT_tuning / max_spike_rate;

CT_tuning = zscore(CT_tuning,0,'all');

diffs = zeros(size(CT_tuning,1)-1,1);
CT0 = CT_tuning(1,:);

[~,I] = max(CT0);
window = I-4:I+4;
window(window<1) = [];
window(window>17) = [];

for ct = 1:size(CT_tuning,1)

    diffs(ct) = trapz(CT0(window)) - trapz(CT_tuning(ct,window));
%     diffs(ct) = mean(CT0(window) - CT_tuning(ct,window));
    
end

p = polyfit(1:4,diffs(2:5),1);

% if p(1) > 0.5
%     sensitivity = 1;
% else
%     sensitivity = 0;
% end

sensitivity = p(1);

end


function nUnits = count_units(unit_list)

nUnits = 0;
for pen = 1:length(unit_list)

    nUnits = nUnits + size(unit_list{pen,3},1);
end

end

% function sensitivity = estimate_pitch_sensitivity(CT_tuning)
% 
% % CT_tuning should be 5 x 17 where each row is the tuning curve for one of
% % the CT stimuli IN ORDER so that CT_tuning(1,:) = tuning for CT0,
% % CT_tuning(2,:) = tuning for CT5, etc up to CT40
% 
% % max_spike_rate = max(max(CT_tuning));
% % CT_tuning = CT_tuning / max_spike_rate;
% 
% CT_tuning = zscore(CT_tuning,0,'all');
% 
% diffs = zeros(size(CT_tuning,1)-1,1);
% CT0 = CT_tuning(1,:);
% 
% [~,I] = max(CT0);
% window = I-4:I+4;
% window(window<1) = [];
% window(window>17) = [];
% 
% for ct = 2:size(CT_tuning,1)
% 
%     diff = trapz(CT_tuning(ct-1,window)) - trapz(CT_tuning(ct,window));
%     diffs(ct) = diff;
% %     diffs(ct) = (diff / (trapz(CT_tuning(ct-1,window))+trapz(CT_tuning(ct,window)))/2) * 100;
%     
% end
% 
% % p = polyfit(1:4,diffs(2:5),1);
% % sensitivity = p(1);
% 
% sensitivity = mean(diffs);
% end