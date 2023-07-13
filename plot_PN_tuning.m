% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, May 2023

load('PNs_05_99')

colors = [1 0 0; 0 0 1; 0 0 0];

figure;
for pen = 1:length(PN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' PN_units{pen,1} '/tmp02/Spikes_' PN_units{pen,1} '_' PN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'high','low','CT0'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    PNUnits = PN_units{pen,3};

    window = [0 0.1];

    for uu = 1:size(PNUnits,1)

        nexttile
        plot_tuning_by_cond(Y,type,F0,PNUnits(uu,1),stims,zeros(length(stims),1),PN_units{pen,1},PN_units{pen,2},window, colors);
        xlabel('')
        ylabel('')
        xticks([])
        yticks([])
    end
end % ends recording loop


sgtitle('')