% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, May 2023

load('HNs_05_99_3.mat')
PN_units = HN_units;
clear HN_units

% colors = [1 0 0; 0 0 1; 0 0 0];
colors = [1 0 0; 1 0 1; 0 0 0; 0 1 1; 0 0 1];

figure;
for pen = 1:length(PN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' PN_units{pen,1} '/final/Spikes_' PN_units{pen,1} '_' PN_units{pen,2} '_Good_Pitch.mat']);

%     stims = {'high','alt','rand'};
    stims = {'CT0','CT5','CT10','CT20','CT40'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    PNUnits = PN_units{pen,3};

    window = [0 0.1];

    for uu = 1:size(PNUnits,1)

        clf
        plot_tuning_by_cond(Y,type,F0,PNUnits(uu,1),stims,zeros(length(stims),1),PN_units{pen,1},PN_units{pen,2},window, colors);
        xlabel('')
        ylabel('')
        xticks([])
        yticks([])

        pause
    end
end % ends recording loop


sgtitle('')