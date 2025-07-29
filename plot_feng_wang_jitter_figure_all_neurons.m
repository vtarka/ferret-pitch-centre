%% Plot distribution of firing rates vs jitter level a la Feng & Wang 2017
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, February 2024

file_name = 'TN_units.mat';

other_type_file = 'HN_units.mat';

mat_struct = load(file_name);
mat_cell = struct2cell(mat_struct);
units_by_rec = mat_cell{1};

mat_struct = load(other_type_file);
mat_cell = struct2cell(mat_struct);
other_units_by_rec = mat_cell{1};

window = [0 0.1];

uCounter = 1;

FRs = zeros(1241,3,17);

for pen = 1:length(units_by_rec)


    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/Batch0324/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

%     units = units_by_rec{pen,3};
    units = unique(Y(:,3));
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    stims = {'CT0','CT5','CT10','CT20','CT40'};
%     stims = {'high','alt','rand'};

    HN_units = units_by_rec{pen,3};
    other_units = other_units_by_rec{pen,3};

    for uu = 1:length(units)

        if responsive(uu)
            unit = units(uu);
    
            if ismember(unit,HN_units) || ismember(unit,other_units)
                continue
            end
    
            unitSpikes = Y(Y(:,3)==unit,:); % get spikes for just this unit
    
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
                    meanSpikes(ff) = mean(spikeRates);
    
                end % ends looping through F0s
    
                FRs(uCounter,ss,:) = meanSpikes;
    
            end % ends loop through stims
    
            uCounter = uCounter + 1;
        end
    end
end

FRs(uCounter:end,:,:) = [];
%% isolate firing rates near CT0's best F0

nUnits = size(FRs,1);

norm_best_FRs = zeros(size(FRs,1),length(stims));

for uu = 1:nUnits

    CT0_tuning = squeeze(FRs(uu,1,:));

    [~,I] = max(CT0_tuning);

    offset = 1;
    window = I-offset:I+offset; % look in the window surrounding the best frequency
    window(window<1) = [];
    window(window>17) = [];
    
    mean_tuning = mean(squeeze(FRs(uu,:,window))','omitnan');

    norm_best_FRs(uu,:) = (mean_tuning - min(mean_tuning))/(max(mean_tuning)-min(mean_tuning));

end

figure; boxchart(norm_best_FRs)

xlabel('temporal jitter (%)')
xticklabels({'0','5','10','20','40'})

ylabel('normalized firing rate')

title('temporal Neurons')

set(gca,'fontsize',20)


%% get significance stats

[p,tbl,stats] = anova1(norm_best_FRs);
[c,m,h] = multcompare(stats);

save('all_neurons_feng_weng_stats','c')
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
