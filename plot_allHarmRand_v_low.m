%% Plots Noah's allHarmRand vs low tuning for harmonicity neurons
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, March 2023

HN_units = TN_units;
%%%%%%%%%%%% WINDOW 1 PLOT %%%%%%%%%%%%%
% for each penetration
figure;
for pen = 1:8
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' HN_units{pen,1} '/tmp02/Spikes_' HN_units{pen,1} '_' HN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'low','allHarm','allHarmRand'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    HNUnits = HN_units{pen,3};

    window = [0 0.06];

    for uu = 1:size(HNUnits,1)

        if length(find(HNUnits==HNUnits(uu,1)))<2 && HNUnits(uu,2)==1

            nexttile
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims,zeros(length(stims),1),HN_units{pen,1},HN_units{pen,2},window);
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
for pen = 1:8
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' HN_units{pen,1} '/tmp02/Spikes_' HN_units{pen,1} '_' HN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'low','allHarm','allHarmRand'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    HNUnits = HN_units{pen,3};

    window = [0.06 0.15];

    for uu = 1:size(HNUnits,1)

        if length(find(HNUnits==HNUnits(uu,1)))<2 && HNUnits(uu,2)==2

            nexttile
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims,zeros(length(stims),1),HN_units{pen,1},HN_units{pen,2},window);
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
for pen = 1:8
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' HN_units{pen,1} '/tmp02/Spikes_' HN_units{pen,1} '_' HN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'low','allHarm','allHarmRand'};
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
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims,zeros(length(stims),1),HN_units{pen,1},HN_units{pen,2},window);
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
for pen = 1:8
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' HN_units{pen,1} '/tmp02/Spikes_' HN_units{pen,1} '_' HN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'low','allHarm','allHarmRand'};
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

            subplot(6,2,sp_counter)
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims,zeros(length(stims),1),HN_units{pen,1},HN_units{pen,2},windows(1,:));
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

            subplot(6,2,sp_counter + 1)
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims,zeros(length(stims),1),HN_units{pen,1},HN_units{pen,2},windows(2,:));
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