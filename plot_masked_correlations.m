%% Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, January 2024

%%

file_name = '0225HNs';

other_type_file = '0225TNs';

mat_struct = load(file_name);
mat_cell = struct2cell(mat_struct);
units_by_rec = mat_cell{1};

mat_struct = load(other_type_file);
mat_cell = struct2cell(mat_struct);
other_units_by_rec = mat_cell{1};

stims_to_plot = {'low','F0MaskLow'};

window = [0 0.1];
hn_corrs = [];
all_corrs = [];

confusion_matrix = zeros(17);
BF0s = [];

for pen = 1:length(units_by_rec)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/Batch0324/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

    HN_units = units_by_rec{pen,3};
    other_units = other_units_by_rec{pen,3};

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    for uu = 1:length(units)

        if responsive(uu)

            unit = units(uu);

            if ~ismember(unit,HN_units)
                continue
            end
    
            unitSpikes = Y(Y(:,3)==unit,:); % get spikes for just this unit
    
            tuning_curves = zeros(2,17);
    
            % go through each stim we want to plot
            for ss = 1:length(stims_to_plot)
        
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
    
                tuning_curves(ss,:) = meanSpikes;
            end
    
            if ismember(unit,HN_units)
                hn_corrs(end+1) = corr(tuning_curves(1,:)',tuning_curves(2,:)');

                [~,I1] = max(tuning_curves(1,:));
                [~,I2] = max(tuning_curves(2,:));
                BF0s = [BF0s; [I1 I2]];

                confusion_matrix(I1,I2) = confusion_matrix(I1,I2) + 1;

            elseif ~ismember(unit,other_units)
                all_corrs(end+1) = corr(tuning_curves(1,:)',tuning_curves(2,:)');
            end

        end
    end
end


% [hn_N,hn_edges] = histcounts(hn_corrs,-1:0.1:1);
% 
% [all_N,all_edges] = histcounts(all_corrs,-1:0.1:1);
% 
% hn_N = hn_N / length(hn_corrs);
% all_N = all_N / length(all_corrs);

% figure;
% hold on
% histogram(hn_corrs,-1:0.1:1,'normalization','probability','facecolor','red','facealpha',0.3,'LineWidth',2)
% histogram(all_corrs,-1:0.1:1,'normalization','probability','facecolor','white','facealpha',0.5,'LineWidth',2)
% 
% set(gca,'fontsize',24)
% 
% legend({'temporal neurons','all other neurons'},'location','northwest')
% xlabel('masked and unmasked high harmonics tuning corr.')
% ylabel('% of neurons')

%% offset histograms

figure;
hold on
[HN_N, HN_edges] = histcounts(hn_corrs,-1:0.1:1,'normalization','probability');
[all_N, all_edges] = histcounts(all_corrs,-1:0.1:1,'normalization','probability');

bar(-0.96:0.1:0.94, HN_N,0.25,'facecolor','blue','linewidth',2)
bar(-0.93:0.1:0.97, all_N,0.25,'facecolor','white','linewidth',2)

set(gca,'fontsize',24)

xticks(-1:0.5:1)

ylim([0 0.25])

legend({'harmonicity neurons','all other neurons'},'location','northwest')
xlabel('masked and unmasked low harmonics tuning corr.')
ylabel('% of neurons')