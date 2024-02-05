%% Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, January 2024

%%

file_name = 'HN_units_new_05';

other_type_file = 'TN_units_new_05';

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

for pen = 1:length(units_by_rec)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/final/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

    HN_units = units_by_rec{pen,3};
    other_units = other_units_by_rec{pen,3};

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    for uu = 1:length(units)

        if responsive(uu)

            unit = units(uu);
    
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

figure;
hold on
histogram(hn_corrs,-1:0.1:1,'normalization','probability','facecolor','black','facealpha',0.3)
histogram(all_corrs,-1:0.1:1,'normalization','probability','facecolor','white','facealpha',0.3)

set(gca,'fontsize',24)

legend({'harmonicity neurons','all other neurons'},'location','northwest')
xlabel('masked and unmasked low harmonics tuning corr.')
ylabel('% of neurons')