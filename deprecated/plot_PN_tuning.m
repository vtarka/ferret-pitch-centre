% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

for ap = 1:length(PN_units)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' PN_units{ap,1} '/tmp/Spikes_' PN_units{ap,1} '_' PN_units{ap,2} '_Good_Pitch.mat']);

    % function a = plot_tuning_by_cond(Y,type,F0,unit,stims,BFs,animal,pen)

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));
    PNUnits = PN_units{ap,3};
    [~,PNUnit_IDXs] = ismember(PNUnits,allUnits);
    BFs = BFs(PNUnit_IDXs,[13,8,1]);

    stims_to_plot = {'tone','allHarm','CT0'};

    for uu = 1:length(PNUnits)
        unit = PNUnits(uu);

        plot_tuning_by_cond(Y,type,F0,unit,stims_to_plot,BFs(uu,:),Animals{ap},Pens{ap});
        pause
    end

end