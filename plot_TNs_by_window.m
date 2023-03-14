%% Produces series of figures showing periodicity (temporal) neurons across 3 response windows
% DEPENDENCIES:
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, March 2023


%%%%%%%%%%%%%%%%%%%%% HIGH LOW CT0 PLOTS %%%%%%%%%%%%%%%%%%

load('TNs_fNoah.mat')

%%
%%%%%%%%%%%% WINDOW 1 PLOT %%%%%%%%%%%%%
% for each penetration
figure;
for pen = 1:length(TN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' TN_units{pen,1} '/tmp02/Spikes_' TN_units{pen,1} '_' TN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'high','low','CT0'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    TNUnits = TN_units{pen,3};

    window = [0 0.06];

    for uu = 1:size(TNUnits,1)

        if length(find(TNUnits==TNUnits(uu,1)))<2 && TNUnits(uu,2)==1

            nexttile
            plot_tuning_by_cond(Y,type,F0,TNUnits(uu,1),stims,zeros(length(stims),1),TN_units{pen,1},TN_units{pen,2},window);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])
        end

    end
end % ends recording loop

sgtitle('Fast onset ONLY sensitive TNs','fontsize',32)

%%
%%%%%%%%%%%% WINDOW 2 PLOT %%%%%%%%%%%%%
% for each penetration
figure;
for pen = 1:length(TN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' TN_units{pen,1} '/tmp02/Spikes_' TN_units{pen,1} '_' TN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'high','low','CT0'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    TNUnits = TN_units{pen,3};

    window = [0.06 0.15];

    for uu = 1:size(TNUnits,1)

        if length(find(TNUnits==TNUnits(uu,1)))<2 && TNUnits(uu,2)==2

            nexttile
            plot_tuning_by_cond(Y,type,F0,TNUnits(uu,1),stims,zeros(length(stims),1),TN_units{pen,1},TN_units{pen,2},window);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])
        end

    end
end % ends recording loop

sgtitle('Slow onset ONLY sensitive TNs','fontsize',32)

%%
%%%%%%%%%%%% WINDOW 3 PLOT %%%%%%%%%%%%%
% for each penetration
figure;
for pen = 1:length(TN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' TN_units{pen,1} '/tmp02/Spikes_' TN_units{pen,1} '_' TN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'high','low','CT0'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    TNUnits = TN_units{pen,3};

    if pen<9
        window = [0.3 0.4];
    else
        window = [0.2 0.3];
    end

    for uu = 1:size(TNUnits,1)

        if length(find(TNUnits==TNUnits(uu,1)))<2 && TNUnits(uu,2)==3

            nexttile
            plot_tuning_by_cond(Y,type,F0,TNUnits(uu,1),stims,zeros(length(stims),1),TN_units{pen,1},TN_units{pen,2},window);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])
        end

    end
end % ends recording loop

sgtitle('Offset ONLY sensitive TNs','fontsize',32)

%%
%%%%%%%%%%%% WINDOW 1 & 2 PLOT %%%%%%%%%%%%%
% for each penetration
figure;
sp_counter = 1;
for pen = 1:length(TN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' TN_units{pen,1} '/tmp02/Spikes_' TN_units{pen,1} '_' TN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'high','low','CT0'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    TNUnits = TN_units{pen,3};

    windows = [0 0.06; 0.06 0.15];

    skip_next = 0;
    for uu = 1:size(TNUnits,1)

        if length(find(TNUnits==TNUnits(uu,1)))==2 && ~skip_next

            if sp_counter > 24
                figure;
                sp_counter = 1;
            end

            subplot(6,4,sp_counter)
            plot_tuning_by_cond(Y,type,F0,TNUnits(uu,1),stims,zeros(length(stims),1),TN_units{pen,1},TN_units{pen,2},windows(1,:));
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

            subplot(6,4,sp_counter + 1)
            plot_tuning_by_cond(Y,type,F0,TNUnits(uu,1),stims,zeros(length(stims),1),TN_units{pen,1},TN_units{pen,2},windows(2,:));
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

sgtitle('')


%%
%%%%%%%%%%%%%%% PHASE SHIFT PLOTS %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WINDOW 1 HIGH ALT RAND PLOT %%%%%%%%%%%%%
% for each penetration
figure;
for pen = 1:length(TN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' TN_units{pen,1} '/tmp02/Spikes_' TN_units{pen,1} '_' TN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'high','alt','rand'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    TNUnits = TN_units{pen,3};

    window = [0 0.06];

    for uu = 1:size(TNUnits,1)

        if length(find(TNUnits==TNUnits(uu,1)))<2 && TNUnits(uu,2)==1

            nexttile
            plot_tuning_by_cond(Y,type,F0,TNUnits(uu,1),stims,zeros(length(stims),1),TN_units{pen,1},TN_units{pen,2},window);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

        end
    end
end

sgtitle('Fast onset ONLY sensitive TNs','fontsize',32)

%%
%%%%%%%%%%%% WINDOW 2 HIGH ALT RAND PLOT %%%%%%%%%%%%%
% for each penetration
figure;
for pen = 1:length(TN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' TN_units{pen,1} '/tmp02/Spikes_' TN_units{pen,1} '_' TN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'high','alt','rand'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    TNUnits = TN_units{pen,3};

    window = [0.06 0.15];

    for uu = 1:size(TNUnits,1)

        if length(find(TNUnits==TNUnits(uu,1)))<2 && TNUnits(uu,2)==2

            nexttile
            plot_tuning_by_cond(Y,type,F0,TNUnits(uu,1),stims,zeros(length(stims),1),TN_units{pen,1},TN_units{pen,2},window);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])
        end
    end
end

sgtitle('Slow onset ONLY sensitive TNs','fontsize',32)

%%
%%%%%%%%%%%% WINDOW 3 HIGH ALT RAND PLOT %%%%%%%%%%%%%
% for each penetration
figure;
for pen = 1:length(TN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' TN_units{pen,1} '/tmp02/Spikes_' TN_units{pen,1} '_' TN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'high','alt','rand'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    TNUnits = TN_units{pen,3};

    if pen<9
        window = [.3 .4];
    else
        window = [.2 .3];
    end

    for uu = 1:size(TNUnits,1)

        if length(find(TNUnits==TNUnits(uu,1)))<2 && TNUnits(uu,2)==3

            nexttile
            plot_tuning_by_cond(Y,type,F0,TNUnits(uu,1),stims,zeros(length(stims),1),TN_units{pen,1},TN_units{pen,2},window);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])
        end
    end
end

sgtitle('Offset ONLY sensitive TNs','fontsize',32)

%%
%%%%%%%%%%%% WINDOW 1 & 2 PLOT %%%%%%%%%%%%%
% for each penetration
figure;
sp_counter = 1;
nPlots = 4;
flag = 1;
for pen = 1:length(TN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' TN_units{pen,1} '/tmp02/Spikes_' TN_units{pen,1} '_' TN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'high','alt','rand'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    TNUnits = TN_units{pen,3};

    windows = [0 0.06; 0.06 0.15];

    
    skip_next = 0;
    for uu = 1:size(TNUnits,1)

        if length(find(TNUnits==TNUnits(uu,1)))==2 && ~skip_next

            if sp_counter > 16 && flag==1
                figure;
                sp_counter = 1;
                nPlots = 5;
                flag = 0;
            end

            subplot(nPlots,4,sp_counter)
            plot_tuning_by_cond(Y,type,F0,TNUnits(uu,1),stims,zeros(length(stims),1),TN_units{pen,1},TN_units{pen,2},windows(1,:));
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

            subplot(nPlots,4,sp_counter + 1)
            plot_tuning_by_cond(Y,type,F0,TNUnits(uu,1),stims,zeros(length(stims),1),TN_units{pen,1},TN_units{pen,2},windows(2,:));
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

sgtitle('')

%%
%%%%%%%%%%%% WINDOW 1 & 2 & 3 PLOT %%%%%%%%%%%%%
% for each penetration
figure;
sp_counter = 1;
nPlots = 7;
for pen = 1:length(TN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' TN_units{pen,1} '/tmp02/Spikes_' TN_units{pen,1} '_' TN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'high','alt','rand'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    TNUnits = TN_units{pen,3};

    if pen < 9
        windows = [0 0.06; 0.06 0.15; 0.3 0.4];
    else
        windows = [0 0.06; 0.06 0.15; 0.2 0.3];
    end

    skip_next = 0;
    for uu = 1:size(TNUnits,1)

        if length(find(TNUnits==TNUnits(uu,1)))==3 && ~skip_next

%             if sp_counter > 24
%                 figure;
%                 sp_counter = 1;
%             end

            subplot(nPlots,3,sp_counter)
            plot_tuning_by_cond(Y,type,F0,TNUnits(uu,1),stims,zeros(length(stims),1),TN_units{pen,1},TN_units{pen,2},windows(1,:));
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

            subplot(nPlots,3,sp_counter + 1)
            plot_tuning_by_cond(Y,type,F0,TNUnits(uu,1),stims,zeros(length(stims),1),TN_units{pen,1},TN_units{pen,2},windows(2,:));
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

            subplot(nPlots,3,sp_counter + 2)
            plot_tuning_by_cond(Y,type,F0,TNUnits(uu,1),stims,zeros(length(stims),1),TN_units{pen,1},TN_units{pen,2},windows(3,:));
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

sgtitle('')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLICK TRAIN PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% for each penetration
figure;
for pen = 1:length(TN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' TN_units{pen,1} '/tmp02/Spikes_' TN_units{pen,1} '_' TN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'CT0','CT5','CT10','CT20','CT40'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    TNUnits = TN_units{pen,3};

    window = [0 0.06];

    for uu = 1:size(TNUnits,1)

        if length(find(TNUnits==TNUnits(uu,1)))<2 && TNUnits(uu,2)==1

            nexttile
            plot_tuning_by_cond(Y,type,F0,TNUnits(uu,1),stims,zeros(length(stims),1),TN_units{pen,1},TN_units{pen,2},window);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

        end
    end
end

sgtitle('Fast onset ONLY sensitive TNs','fontsize',32)

%%
%%%%%%%%%%%% WINDOW 2 CT PLOT %%%%%%%%%%%%%
% for each penetration
figure;
for pen = 1:length(TN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' TN_units{pen,1} '/tmp02/Spikes_' TN_units{pen,1} '_' TN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'CT0','CT5','CT10','CT20','CT40'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    TNUnits = TN_units{pen,3};

    window = [0.06 0.15];

    for uu = 1:size(TNUnits,1)

        if length(find(TNUnits==TNUnits(uu,1)))<2 && TNUnits(uu,2)==2

            nexttile
            plot_tuning_by_cond(Y,type,F0,TNUnits(uu,1),stims,zeros(length(stims),1),TN_units{pen,1},TN_units{pen,2},window);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])
        end
    end
end

sgtitle('Slow onset ONLY sensitive TNs','fontsize',32)


%%
%%%%%%%%%%%% WINDOW 3 CT PLOT %%%%%%%%%%%%%
% for each penetration
figure;
for pen = 1:length(TN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' TN_units{pen,1} '/tmp02/Spikes_' TN_units{pen,1} '_' TN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'CT0','CT5','CT10','CT20','CT40'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    TNUnits = TN_units{pen,3};

    if pen<9
        window = [.3 .4];
    else
        window = [.2 .3];
    end

    for uu = 1:size(TNUnits,1)

        if length(find(TNUnits==TNUnits(uu,1)))<2 && TNUnits(uu,2)==3

            nexttile
            plot_tuning_by_cond(Y,type,F0,TNUnits(uu,1),stims,zeros(length(stims),1),TN_units{pen,1},TN_units{pen,2},window);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])
        end
    end
end

sgtitle('Offset ONLY sensitive TNs','fontsize',32)


%%
%%%%%%%%%%%% WINDOW 1 & 2 PLOT %%%%%%%%%%%%%
% for each penetration
figure;
sp_counter = 1;
nPlots = 4;
flag = 1;
for pen = 1:length(TN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' TN_units{pen,1} '/tmp02/Spikes_' TN_units{pen,1} '_' TN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'CT0','CT5','CT10','CT20','CT40'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    TNUnits = TN_units{pen,3};

    windows = [0 0.06; 0.06 0.15];

    
    skip_next = 0;
    for uu = 1:size(TNUnits,1)
        
        if length(find(TNUnits==TNUnits(uu,1)))==2 && ~skip_next

            if sp_counter > 16 && flag==1
                figure;
                sp_counter = 1;
                nPlots = 5;
                flag = 0;
            end

            subplot(nPlots,4,sp_counter)
            plot_tuning_by_cond(Y,type,F0,TNUnits(uu,1),stims,zeros(length(stims),1),TN_units{pen,1},TN_units{pen,2},windows(1,:));
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

            subplot(nPlots,4,sp_counter + 1)
            plot_tuning_by_cond(Y,type,F0,TNUnits(uu,1),stims,zeros(length(stims),1),TN_units{pen,1},TN_units{pen,2},windows(2,:));
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

sgtitle('')

%%
%%%%%%%%%%%% WINDOW 1 & 2 & 3 PLOT %%%%%%%%%%%%%
% for each penetration
figure;
sp_counter = 1;
nPlots = 7;
for pen = 1:length(TN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' TN_units{pen,1} '/tmp02/Spikes_' TN_units{pen,1} '_' TN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'CT0','CT5','CT10','CT20','CT40'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    TNUnits = TN_units{pen,3};

    if pen < 9
        windows = [0 0.06; 0.06 0.15; 0.3 0.4];
    else
        windows = [0 0.06; 0.06 0.15; 0.2 0.3];
    end

    skip_next = 0;
    for uu = 1:size(TNUnits,1)

        if length(find(TNUnits==TNUnits(uu,1)))==3 && ~skip_next

%             if sp_counter > 24
%                 figure;
%                 sp_counter = 1;
%             end

            subplot(nPlots,3,sp_counter)
            plot_tuning_by_cond(Y,type,F0,TNUnits(uu,1),stims,zeros(length(stims),1),TN_units{pen,1},TN_units{pen,2},windows(1,:));
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

            subplot(nPlots,3,sp_counter + 1)
            plot_tuning_by_cond(Y,type,F0,TNUnits(uu,1),stims,zeros(length(stims),1),TN_units{pen,1},TN_units{pen,2},windows(2,:));
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

            subplot(nPlots,3,sp_counter + 2)
            plot_tuning_by_cond(Y,type,F0,TNUnits(uu,1),stims,zeros(length(stims),1),TN_units{pen,1},TN_units{pen,2},windows(3,:));
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

sgtitle('')


%% Plot TN_Unit locations

loc_frequencies = [];
loc_total_units = zeros(5,1);

pen_frequencies = [];
pen_total_units = zeros(20,1);

locs = [1 1 2 1 5 3 5 4 3 3 1 2 3 2 1 1 1 3 4 3];

for pen = 1:length(TN_units)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' TN_units{pen,1} '/tmp02/Spikes_' TN_units{pen,1} '_' TN_units{pen,2} '_Good_Pitch.mat']);
    
    TNUs = TN_units{pen,3};

    loc_total_units(locs(pen)) = loc_total_units(locs(pen)) + length(unique(Y(:,3)));
    pen_total_units(pen) = length(unique(Y(:,3)));

    for uu = 1:size(TNUs,1)

        loc_frequencies = [loc_frequencies; TNUs(uu,2) locs(pen)];
        pen_frequencies = [pen_frequencies; TNUs(uu,2) pen];

    end
end

w1_idx = loc_frequencies(:,1)==1;
w2_idx = loc_frequencies(:,1)==2;
w3_idx = loc_frequencies(:,1)==3;

fsz = 22;

figure;

%%%%%%%%%%%%%%%%%% FAST ONSET %%%%%%%%%%%%%%%%%%%%
subplot(2,3,1)
w1_pens = zeros(4,7);
c = histcounts(pen_frequencies(w1_idx,2),1:21);

for i = 1:4
    w1_pens(i,1:length(find(locs==i))) = sort(c(locs==i),'descend');
end

H = bar(1:4,w1_pens,'stacked');

xticks(1:4)
xticklabels({'lA1','hA1','lAAF','hAAF'})
ylim([0 15])
ylabel('# of units')
set(gca,'fontsize',fsz)
title('Fast Onset')

subplot(2,3,4)
w1_pens_norm = zeros(4,7);
nc = c ./ pen_total_units';

for i = 1:4
    w1_pens_norm(i,1:length(find(locs==i))) = sort(nc(locs==i),'descend');
end

H_norm = bar(1:4,w1_pens_norm,'stacked');

xticks(1:4)
xticklabels({'lA1','hA1','lAAF','hAAF'})
ylim([0 0.25])
ylabel('# of units / total units in pen')
set(gca,'fontsize',fsz)
    

%%%%%%%%%%%%%%%%%% SLOW ONSET %%%%%%%%%%%%%%%%%%%%
subplot(2,3,2)
w2_pens = zeros(4,7);
c = histcounts(pen_frequencies(w2_idx,2),1:21);

for i = 1:4
    w2_pens(i,1:length(find(locs==i))) = sort(c(locs==i),'descend');
end

H = bar(1:4,w2_pens,'stacked');

xticks(1:4)
xticklabels({'lA1','hA1','lAAF','hAAF'})
ylim([0 15])
ylabel('# of units')
set(gca,'fontsize',fsz)
title('Slow Onset')

subplot(2,3,5)
w2_pens_norm = zeros(4,7);
nc = c ./ pen_total_units';

for i = 1:4
    w2_pens_norm(i,1:length(find(locs==i))) = sort(nc(locs==i),'descend');
end

H_norm = bar(1:4,w2_pens_norm,'stacked');

xticks(1:4)
xticklabels({'lA1','hA1','lAAF','hAAF'})
ylim([0 0.25])
ylabel('# of units / total units in pen')
set(gca,'fontsize',fsz)


%%%%%%%%%%%%%%%%%% OFFSET %%%%%%%%%%%%%%%%%%%%
subplot(2,3,3)
w3_pens = zeros(4,7);
c = histcounts(pen_frequencies(w3_idx,2),1:21);

for i = 1:4
    w3_pens(i,1:length(find(locs==i))) = sort(c(locs==i),'descend');
end

H = bar(1:4,w3_pens,'stacked');

xticks(1:4)
xticklabels({'lA1','hA1','lAAF','hAAF'})
ylim([0 15])
ylabel('# of units')
set(gca,'fontsize',fsz)
title('Offset')

subplot(2,3,6)
w3_pens_norm = zeros(4,7);
nc = c ./ pen_total_units';

for i = 1:4
    w3_pens_norm(i,1:length(find(locs==i))) = sort(nc(locs==i),'descend');
end

H_norm = bar(1:4,w3_pens_norm,'stacked');

xticks(1:4)
xticklabels({'lA1','hA1','lAAF','hAAF'})
ylim([0 0.25])
ylabel('# of units / total units in pen')
set(gca,'fontsize',fsz)