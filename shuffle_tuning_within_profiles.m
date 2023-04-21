%% Test how well-aligned the tuning is across stimuli by shuffling responses within a tuning curve
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, March 2023


load('HNs_fNoah')
units_by_rec = HN_units;

stims = {'CT0','CT5','CT10','allHarm','low'};

nNullRuns = 1000;

nUnits = count_units(units_by_rec);
corrs = zeros(nUnits,2);
unit_counter = 1;
real_profiles = zeros(nUnits,length(stims),17);

locs = [1 1 2 1 5 3 5 4 3 3 1 2 3 2 1 1 1 3 4 3];
loc_list = zeros(nUnits,1);

zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));

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
            profile(ss,:) = zscore(meanSpikes);

            if strcmp('low',stims{ss})
                profile(ss,1:2) = nan;
            end

        end % ends stim loop

        corrs(unit_counter,1) = get_avg_pairwise_corr(profile);
        real_profiles(unit_counter,:,:) = profile;

        shuffled_profile = profile;
        shuffled_corrs = zeros(nNullRuns,1);
        for i = 1:nNullRuns
            for ss = 1:length(stims)
                
                shuffle_key = randperm(size(profile,2)); % random permutation of the number of frequencies we presented

                shuffled_profile(ss,:) = profile(ss,shuffle_key);
            end

            shuffled_corrs(i) = get_avg_pairwise_corr(shuffled_profile);

        end

        null_grand_avg_corr = prctile(shuffled_corrs,99.9);
        corrs(unit_counter,2) = null_grand_avg_corr;
        loc_list(unit_counter) = locs(pen);

        unit_counter = unit_counter + 1;
    end
end

%%

% all_corrs = [null_grand_avg_corrs; real_corrs];
% g1 = repmat({'Null'},length(null_grand_avg_corrs),1);
% g2 = repmat({'Real'},length(real_corrs),1);
% g = [g1;g2];
% figure; boxplot(all_corrs,g)
% ylabel('Correlation')
% set(gca,'fontsize',24)
% yline(median(real_corrs),'r--','linewidth',2)
% yline(mean(real_corrs),'k','linewidth',2)


passing_units = find(corrs(:,1)>corrs(:,2));
passing_locations = [];

for i = 1:length(corrs)
    if ismember(i,passing_units)

        [U,S,V] = svd(squeeze(real_profiles(i,:,:)));
        s2_reconstruction = U(:,1:2)*S(1:2,1:2)*V(:,1:2)';
        figure(1); nexttile;

        imagesc(s2_reconstruction)
        axis('off')
% 
%         imagesc(squeeze(real_profiles(i,:,:)));
%         xticks([])
%         yticks([])
    else

        [U,S,V] = svd(squeeze(real_profiles(i,:,:)));
        s2_reconstruction = U(:,1:2)*S(1:2,1:2)*V(:,1:2)';
        figure(2); nexttile;
        imagesc(s2_reconstruction)
        axis('off')

%         imagesc(squeeze(real_profiles(i,:,:)));
%         xticks([])
%         yticks([])
    end
end

figure(1); sgtitle('Passing Pitch Units','fontsize',28)
figure(2); sgtitle('Discarded Pitch Units','fontsize',28)

ymax = 19;
figure;
subplot(1,2,1); histogram(loc_list(passing_units))
xlim([0 4.5])
xticks(1:4)
xticklabels({'lA1','hA1','lAAF','hAAF'})
ylim([0 ymax])
yticks(0:2:ymax)
title('Passing HNs')
set(gca,'fontsize',28)

subplot(1,2,2); histogram(loc_list(corrs(:,1)<corrs(:,2)));
xlim([0 4.5])
xticks(1:4)
xticklabels({'lA1','hA1','lAAF','hAAF'})
ylim([0 ymax])
yticks(0:2:ymax)
title('Not Passing HNs')
set(gca,'fontsize',28)

function c = get_avg_pairwise_corr(tuning)

    % 1st dimension of tuning is the different stimuli
    % 2nd dimension of tuning is the z-scored response across the 17 frequencies

    nStims = size(tuning,1);
    corrs = zeros(nchoosek(nStims,2),1);

    corr_counter = 1;
    for i = 1:nStims
        for k = 1:nStims

            if i<k
                corrs(corr_counter) = corr(tuning(i,:)',tuning(k,:)','rows','complete');
                corr_counter = corr_counter + 1;
            end
        end
    end

    c = mean(corrs,'omitnan');
end


function nUnits = count_units(unit_list)

nUnits = 0;
for pen = 1:length(unit_list)

    nUnits = nUnits + size(unit_list{pen,3},1);
end

end