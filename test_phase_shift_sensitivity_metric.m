%% Plot phase shift tuning ordered by phase shift sensitivity to test given metric
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, April 2023

load('TNs_05_99')
units_by_rec = TN_units;

stims = {'high','alt','rand','F0MaskHigh','tone'};

nUnits = count_units(units_by_rec);
unit_counter = 1;
phase_profiles = zeros(nUnits,length(stims),17);
phase_sensitivity = zeros(nUnits,1);
locs = [1 1 2 1 5 3 5 4 3 3 1 2 3 2 1 1 1 3 4 3];
loc_list = zeros(nUnits,1);

for pen = 1:length(units_by_rec)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/tmp02/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

    units = units_by_rec{pen,3};
% 
%     if ~isempty(units)
%         units(units(:,2)==3,:) = [];
%         units = unique(units(:,1));
%     end

    Flist = unique(F0);
    repeats = unique(Y(:,5));

    for uu = 1:size(units,1)

        unit = units(uu,1);

        unitSpikes = Y(Y(:,3)==unit,:); % get spikes for just this unit
        profile = zeros(length(stims),17);

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

        window = [0 0.1];
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
            profile(ss,:) = (meanSpikes);

            if strcmp('low',stims{ss})
                profile(ss,1:2) = nan;
            end

        end % ends stim loop

        phase_sensitivity(unit_counter) = estimate_phase_shift_sensitivity(profile(1:3,:));
        phase_profiles(unit_counter,:,:) = profile;
        loc_list(unit_counter) = locs(pen);

        unit_counter = unit_counter + 1;
    end
end

[B,I] = sort(phase_sensitivity);
c3 = [0 246 148; 34 0 255; 246 164 0; 0 0 0]/255;

figure;
for uu = 1:nUnits
    nexttile; hold on;
    for ss = 1:length(stims)-2
        plot(1:17,squeeze(phase_profiles(I(uu),ss,:)),'Color',c3(ss,:),'linewidth',2)
    end
    axis tight
    xticks([])
    yticks([])
    if B(uu)==1
        title('Sensitive')
    else
        title('Insensitive')
    end
    
end

%%
figure;
specialUs = [6 8 11 15 48 51 52];
for uu = 1:length(specialUs)

    subplot(length(specialUs),3,(uu*3)-2); hold on
    for ss = 1:3
        plot(1:17,squeeze(phase_profiles(specialUs(uu),ss,:)),'Color',c3(ss,:),'linewidth',2)
    end
    axis tight

    axis tight
    xticks([])
    yticks([])


    subplot(length(specialUs),3,(uu*3)-1); hold on
    special_order = [4 2 3];
    for ss = 1:3
        plot(1:17,squeeze(phase_profiles(specialUs(uu),special_order(ss),:)),'Color',c3(ss,:),'linewidth',2)
    end
    axis tight

    axis tight
    xticks([])
    yticks([])

    subplot(length(specialUs),3,(uu*3)); hold on
    special_order = [1 2 3 5];
    for ss = 1:4
        if ss == 4
            plot(1:17,squeeze(phase_profiles(specialUs(uu),special_order(ss),:)),'Color',c3(ss,:),'linewidth',3,'linestyle','--')
        else
            plot(1:17,squeeze(phase_profiles(specialUs(uu),special_order(ss),:)),'Color',c3(ss,:),'linewidth',2)
        end
    end
    axis tight

    axis tight
    xticks([])
    yticks([])
    
end

sgtitle('Phase Sensitive and Insensitive Temporal Neurons','fontsize',28)

ymax = 12;
figure; subplot(1,2,1)
histogram(loc_list(phase_sensitivity==0))
xticks(1:4)
xticklabels({'lA1','hA1','lAAF','hAAF'})
ylim([0 ymax])
yticks(0:2:ymax)
title('Not Phase Sensitive TNs')
set(gca,'fontsize',24)

subplot(1,2,2)
histogram(loc_list(phase_sensitivity==1))
xticks(1:4)
xticklabels({'lA1','hA1','lAAF','hAAF'})
ylim([0 ymax])
yticks(0:2:ymax)
title('Phase Sensitive TNs')
set(gca,'fontsize',24)


% 
% figure;
% n_subplots = length(non_masked_sensitive*2);
% 
% for uu = 1:length(non_masked_sensitive)
%     subplot(length(non_masked_sensitive),2,(uu*2)-1)
%     hold on
%     for ss = 1:length(stims)
%         plot(1:17,squeeze(non_masked_profiles(non_masked_sensitive(uu),ss,:)),'Color',colors(ss,:))
%     end
%     axis tight
% 
%     subplot(length(non_masked_sensitive),2,uu*2)
%     hold on
%     for ss = 1:length(stims)
%         plot(1:17,squeeze(phase_profiles(non_masked_sensitive(uu),ss,:)),'Color',colors(ss,:))
%     end
%     axis tight
% end


% function sensitivity = estimate_phase_shift_sensitivity(phase_tuning)
% 
% [r_alt,lags_alt] = xcorr(phase_tuning(1,:),phase_tuning(2,:));
% [~,max_r_idx] = max(r_alt);
% best_lag_alt = lags_alt(max_r_idx);
% 
% [r_rand,lags_rand] = xcorr(phase_tuning(1,:),phase_tuning(3,:));
% [~,max_r_idx] = max(r_rand);
% best_lag_rand = lags_rand(max_r_idx);
% 
% zscored_tuning = zscore(phase_tuning,0,'all');
% rand_zscored_tuning = zscored_tuning(3,:);
% 
% 
% if best_lag_alt > 0 && (max(rand_zscored_tuning) < 0.5 || best_lag_rand > 0)
%     sensitivity = 1;
% else
%     sensitivity = 0;
% end
% 
% 
% % if max(rand_zscored_tuning)<0.5 
% %     if best_lag_alt>0
% %         sensitivity = 1;
% %     else
% %         sensitivity = 0;
% %     end
% % else 
% %     if best_lag_rand > 0
% %         if best_lag_alt > 0
% %             sensitivity = 1;
% %         else
% %             sensitivity = -1;
% %         end
% %     else
% %         sensitivity = -1;
% %     end
% % end
% 
% end


% function nUnits = count_units(unit_list)
% 
% nUnits = 0;
% for pen = 1:length(unit_list)
% 
%     nUnits = nUnits + size(unit_list{pen,3},1);
% end
% 
% end


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