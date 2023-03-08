%% Produces series of figures showing harmonicity neurons across 3 response windows
% DEPENDENCIES:
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, March 2023

%%%%%%%%%%%% WINDOW 1 PLOT %%%%%%%%%%%%%
% for each penetration
figure;
for pen = 1:length(HN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' HN_units{pen,1} '/tmp/Spikes_' HN_units{pen,1} '_' HN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'high','low','CT0'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    HNUnits = HN_units{pen,3};

    window = [0 0.06];

    for uu = 1:size(HNUnits,1)

        if length(find(HNUnits==HNUnits(uu,1)))<2 && HNUnits(uu,2)==1

            nexttile
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims,zeros(length(stims),1),Animals{ap},Pens{ap},window);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])
        end

    end
end % ends recording loop

sgtitle('Fast onset ONLY sensitive HNs','fontsize',32)


%%
%%%%%%%%%%%% WINDOW 2 PLOT %%%%%%%%%%%%%
% for each penetration
figure;
for pen = 1:length(HN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' HN_units{pen,1} '/tmp/Spikes_' HN_units{pen,1} '_' HN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'high','low','CT0'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    HNUnits = HN_units{pen,3};

    window = [0.06 0.15];

    for uu = 1:size(HNUnits,1)

        if length(find(HNUnits==HNUnits(uu,1)))<2 && HNUnits(uu,2)==2

            nexttile
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims,zeros(length(stims),1),Animals{ap},Pens{ap},window);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])
        end

    end
end % ends recording loop

sgtitle('Slow onset ONLY sensitive HNs','fontsize',32)

%%
%%%%%%%%%%%% WINDOW 3 PLOT %%%%%%%%%%%%%
% for each penetration
figure;
for pen = 1:length(HN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' HN_units{pen,1} '/tmp/Spikes_' HN_units{pen,1} '_' HN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'high','low','CT0'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    HNUnits = HN_units{pen,3};

    if pen<9
        window = [0.3 0.4];
    else
        window = [0.2 0.3];
    end

    for uu = 1:size(HNUnits,1)

        if length(find(HNUnits==HNUnits(uu,1)))<2 && HNUnits(uu,2)==3

            nexttile
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims,zeros(length(stims),1),Animals{ap},Pens{ap},window);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])
        end

    end
end % ends recording loop

sgtitle('Offset ONLY sensitive HNs','fontsize',32)


%%
%%%%%%%%%%%% WINDOW 1 & 2 PLOT %%%%%%%%%%%%%
% for each penetration
figure;
sp_counter = 1;
for pen = 1:length(HN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' HN_units{pen,1} '/tmp/Spikes_' HN_units{pen,1} '_' HN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'high','low','CT0'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    HNUnits = HN_units{pen,3};

    windows = [0 0.06; 0.06 0.15];

    skip_next = 0;
    for uu = 1:size(HNUnits,1)

        if length(find(HNUnits==HNUnits(uu,1)))==2 && ~skip_next

            if sp_counter > 24
                figure;
                sp_counter = 1;
            end

            subplot(6,4,sp_counter)
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims,zeros(length(stims),1),Animals{ap},Pens{ap},windows(1,:));
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

            subplot(6,4,sp_counter + 1)
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims,zeros(length(stims),1),Animals{ap},Pens{ap},windows(2,:));
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

            sp_counter = sp_counter + 2;

            skip_next = 1;

        else
            skip_next = 0;
        end
    end
end % ends recording loop

% sgtitle('Offset ONLY sensitive HNs','fontsize',32)

%%
%%%%%%%%%%%% WINDOW 1 & 2 & 3 PLOT %%%%%%%%%%%%%
% for each penetration
figure;
sp_counter = 1;
for pen = 1:length(HN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' HN_units{pen,1} '/tmp/Spikes_' HN_units{pen,1} '_' HN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'high','low','CT0'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    HNUnits = HN_units{pen,3};

    if pen < 9
        windows = [0 0.06; 0.06 0.15; 0.3 0.4];
    else
        windows = [0 0.06; 0.06 0.15; 0.2 0.3];
    end

    skip_next = 0;
    for uu = 1:size(HNUnits,1)

        if length(find(HNUnits==HNUnits(uu,1)))==3 && ~skip_next

            if sp_counter > 24
                figure;
                sp_counter = 1;
            end

            subplot(2,3,sp_counter)
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims,zeros(length(stims),1),Animals{ap},Pens{ap},windows(1,:));
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

            subplot(2,3,sp_counter + 1)
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims,zeros(length(stims),1),Animals{ap},Pens{ap},windows(2,:));
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

            subplot(2,3,sp_counter + 2)
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims,zeros(length(stims),1),Animals{ap},Pens{ap},windows(3,:));
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

            sp_counter = sp_counter + 3;

            skip_next = 1;

        else
            skip_next = 0;
        end
    end
end % ends recording loop


%% Plot HN_Unit locations

loc_frequencies = [];

locs = [1 1 2 1 5 3 5 4 3 3 1 2 3 2 1 1 1 3 4 3]; 

for pen = 1:length(HN_units)
    
    HNUs = HN_units{pen,3};

    for uu = 1:size(HNUs,1)

        loc_frequencies = [loc_frequencies; HNUs(uu,2) locs(pen)];
    end

end

w1_idx = loc_frequencies(:,1)==1;
w2_idx = loc_frequencies(:,1)==2;
w3_idx = loc_frequencies(:,1)==3;

figure;
subplot(1,3,1)
histogram(loc_frequencies(w1_idx,2))
xticks(1:4)
xticklabels({'lA1','hA1','lAAF','hAAF'})
ylim([0 25])
set(gca,'fontsize',18)
title('Fast Onset')

subplot(1,3,2)
histogram(loc_frequencies(w2_idx,2))
xlim([0 4.7])
xticks(1:4)
xticklabels({'lA1','hA1','lAAF','hAAF'})
ylim([0 25])
set(gca,'fontsize',18)
title('Slow Onset')

subplot(1,3,3)
histogram(loc_frequencies(w3_idx,2))
xlim([0 4.7])
xticks(1:4)
xticklabels({'lA1','hA1','lAAF','hAAF'})
ylim([0 25])
set(gca,'fontsize',18)
title('Offset')
