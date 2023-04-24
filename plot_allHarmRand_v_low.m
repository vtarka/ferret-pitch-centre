%% Plots Noah's allHarmRand vs low tuning for harmonicity neurons
% DEPENDENCIES: plot_tuning_by_cond.m 
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, March 2023

load('HNs.mat')

% for each penetration
figure;
for pen = 1:8
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' HN_units{pen,1} '/tmp02/Spikes_' HN_units{pen,1} '_' HN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'low','allHarm','allHarmRand'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    HNUnits = HN_units{pen,3};

    window = [0 0.1];

    for uu = 1:length(HNUnits)
        nexttile
        plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims,zeros(length(stims),1),HN_units{pen,1},HN_units{pen,2},window);
        xlabel('')
        ylabel('')
        xticks([])
        yticks([])
    end
end % ends recording loop