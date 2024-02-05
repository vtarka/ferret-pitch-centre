
window = [0 0.1];


%%%%%%%%%%%%%%%%%%% HNs %%%%%%%%%%%%%%%%%%%%%
figure;
load('HN_units_f')
units_by_rec = HN_units;
stims_for_profile = {'CT0','CT5','allHarm','low'};

uCounter = 1;

for pen = 1:length(HN_units)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/tmp02/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

    units = units_by_rec{pen,3};

    for uu = 1:length(units)

        if uCounter == 20 || uCounter == 17
            nexttile;
            imagesc(get_response_profile(Y,type,F0,units(uu),stims_for_profile,window))
            yticks([])
            xticks([])
            title('Lost')
        end

        uCounter = uCounter + 1;
    end
end

load('HN_units_f2')
units_by_rec = HN_units;
stims_for_profile = {'CT0','CT5','allHarm','low'};

uCounter = 1;

for pen = 1:length(HN_units)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/tmp02/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

    units = units_by_rec{pen,3};

    for uu = 1:length(units)

        if uCounter == 28 || uCounter == 25
            nexttile;
            imagesc(get_response_profile(Y,type,F0,units(uu),stims_for_profile,window))
            yticks([])
            xticks([])
            title('Gained')

            
        end

        uCounter = uCounter + 1;
    end
end

sgtitle('Harmonicity Neurons','fontsize',22)

%%%%%%%%%%%%%%%%%%% TNs %%%%%%%%%%%%%%%%%%%%%
figure;
load('TN_units_f')
units_by_rec = TN_units;
stims_for_profile = {'CT0','CT5','allHarm','high'};

uCounter = 1;

for pen = 1:length(TN_units)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/tmp02/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

    units = units_by_rec{pen,3};

    for uu = 1:length(units)

        if uCounter == 34 || uCounter == 25 || uCounter == 10
            nexttile;
            imagesc(get_response_profile(Y,type,F0,units(uu),stims_for_profile,window))
            yticks([])
            xticks([])
            title('Lost')
        end

        uCounter = uCounter + 1;
    end
end

load('TN_units_f2')
units_by_rec = TN_units;
stims_for_profile = {'CT0','CT5','allHarm','high'};

uCounter = 1;

for pen = 1:length(TN_units)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/tmp02/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

    units = units_by_rec{pen,3};

    for uu = 1:length(units)

        if uCounter == 33 || uCounter == 2
            nexttile;
            imagesc(get_response_profile(Y,type,F0,units(uu),stims_for_profile,window))
            yticks([])
            xticks([])
            title('Gained')
        end

        uCounter = uCounter + 1;
    end
end

sgtitle('Temporal Neurons','fontsize',22)



%%%%%%%%%%%%%%%%%%% PNs %%%%%%%%%%%%%%%%%%%%%
figure;
load('PN_units_f2')
units_by_rec = PN_units;
stims_for_profile = {'CT0','CT5','allHarm','high','low'};

uCounter = 1;

for pen = 1:length(TN_units)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/tmp02/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

    units = units_by_rec{pen,3};

    for uu = 1:length(units)

        if uCounter == 31
            nexttile;
            imagesc(get_response_profile(Y,type,F0,units(uu),stims_for_profile,window))
            yticks([])
            xticks([])
            title('Gained')

            figure;
            plot_tuning_by_cond(Y,type,F0,units(uu),{'CT0','F0MaskLow','F0MaskHigh'},[0 0 0],'hi','9',[0 0.1],[0 0 0; 0 0 1; 1 0 0])
        end

        uCounter = uCounter + 1;
    end
end


sgtitle('Pitch Neurons','fontsize',22)
