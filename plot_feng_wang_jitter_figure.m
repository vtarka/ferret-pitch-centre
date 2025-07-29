%% Plot distribution of firing rates vs jitter level a la Feng & Wang 2017
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, February 2024

file_name = '0225TNs.mat';

mat_struct = load(file_name);
mat_cell = struct2cell(mat_struct);
units_by_rec = mat_cell{1};

window = [0 0.1];

nUnits = count_units(units_by_rec);

FRs = zeros(nUnits,3,17);
uCounter = 1;

for pen = 1:length(units_by_rec)


    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/Batch0225/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

    units = units_by_rec{pen,3};
%     units = unique(Y(:,3));
    Flist = unique(F0);
    repeats = unique(Y(:,5));
%     stims = {'CT0','CT5','CT10','CT20','CT40'};
    stims = {'high','alt','rand'};

    for uu = 1:length(units)

        unit = units(uu);

        unitSpikes = Y(Y(:,3)==unit,:); % get spikes for just this unit

        maxRate = 0;

        % for each stim type
        for ss = 1:length(stims)
                
            meanSpikes = zeros(length(Flist),1);

            % go through each F0 to find BF
            for ff = 1:length(Flist)

                stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); % unique name for combination of stim type and F0

                if isempty(stimNum) % if this stim type and fo combo wasn't presented
                    continue
                end

                nSpikes = zeros(length(repeats),1); % initialize space to count the number of spikes per trial repeat

                % for each repeat
                for rr = 1:length(repeats)

                    spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                    nSpikes(rr) = sum(spikeIDXs);

                end % ends repeat loop

                spikeRates = nSpikes / diff(window);

                maxRateContender = max(spikeRates);

                if maxRateContender > maxRate
                    maxRate = maxRateContender;
                end

                meanSpikes(ff) = mean(spikeRates);

            end % ends looping through F0s

            FRs(uCounter,ss,:) = meanSpikes;

        end % ends loop through stims

%         maxRate = max(max(squeeze(FRs(uCounter,:,:))));

        FRs(uCounter,:,:) = FRs(uCounter,:,:);

        uCounter = uCounter + 1;
    end
end

%% isolate firing rates near CT0's best F0

norm_best_FRs = zeros(nUnits,length(stims));

BF0s = zeros(nUnits,5);

for uu = 1:nUnits

    CT0_tuning = squeeze(FRs(uu,1,:));
%     CT5_tuning = squeeze(FRs(uu,2,:));
%     CT10_tuning = squeeze(FRs(uu,3,:));
%     CT20_tuning = squeeze(FRs(uu,4,:));
%     CT40_tuning = squeeze(FRs(uu,5,:));

    [maxR,I] = max(CT0_tuning);
    minR = min(CT0_tuning);

%     [~,I2] = max(CT5_tuning);
%     [~,I3] = max(CT10_tuning);
%     [~,I4] = max(CT20_tuning);
%     [~,I5] = max(CT40_tuning);
% 
%     BF0s(uu,:) = [I I2 I3 I4 I5];

    offset = 2;
    window = I-offset:I+offset; % look in the window surrounding the best frequency
    window(window<1) = [];
    window(window>17) = [];
    
    mean_tuning = median(squeeze(FRs(uu,:,window))','omitnan');
%     mean_tuning = squeeze(FRs(uu,:,window));

    if isnan(mean_tuning(1))
        a=1;
    end
%     norm_best_FRs(uu,:) = (mean_tuning - min(mean_tuning))/(max(mean_tuning)-min(mean_tuning));
    norm_best_FRs(uu,:) = (mean_tuning / mean(mean_tuning));
%     norm_best_FRs(uu,:) = mean_tuning; %(mean_tuning - min(mean_tuning))/(max(mean_tuning)-min(mean_tuning));

end

figure; boxplot(norm_best_FRs);

xlabel('temporal jitter (%)')
xticklabels({'0','5','10','20','40'})

ylabel('normalized firing rate')

title('temporal Neurons')

set(gca,'fontsize',20)


%% get significance stats

[p,tbl,stats] = anova1(norm_best_FRs);
[c,m,h] = multcompare(stats);

% save('HN_feng_weng_stats_phase_2','c')
%%
% norm_FRs = normalize(FRs);

% norm_FRs = (FRs - min(min(FRs)))/(max(max(FRs))-min(min(FRs)));

% for i = 1:nUnits
% 
%     this_jlevel = FRs(i,:);
% 
%     this_jlevel_norm = (this_jlevel - min(this_jlevel))/(max(this_jlevel) - min(this_jlevel));
% 
%     norm_FRs(i,:) = this_jlevel_norm;
% 
% end

% figure; boxchart(norm_FRs);


%% plot best F0 distribution
BF0s_diffs = BF0s;

for uu = 1:nUnits
    BF0s_diffs(uu,2:5) = BF0s(uu,2:5) - BF0s(uu,1);
end

figure; histogram(abs(BF0s_diffs(:,2:5)))
