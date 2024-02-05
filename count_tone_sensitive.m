% Counts how many TNs, HNs, and PNs, are sensitive to pure tones
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, May 2023

load('HNs_05_99.mat')
units_by_rec = HN_units;

tone_sensitive = 0;

for pen = 1:length(units_by_rec)
            
        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/tmp02/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

        stims = unique(type);
        tone_number = find(strcmp(stims,'tone'));
        special_units = units_by_rec{pen,3};
        units = unique(Y(:,3));
        Flist = unique(F0);
        repeats = unique(Y(:,5));
        

        for uu = 1:length(units)

            unit = units(uu);

            if ismember(unit,special_units) && sensitivity(uu,tone_number)
                tone_sensitive = tone_sensitive + 1;
            end
        end
end


%%

load('TNs_05_99.mat')
units_by_rec = TN_units;

tone_sensitive = 0;

for pen = 1:length(units_by_rec)
            
        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/tmp02/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

        stims = unique(type);
        tone_number = find(strcmp(stims,'tone'));
        special_units = units_by_rec{pen,3};
        units = unique(Y(:,3));
        Flist = unique(F0);
        repeats = unique(Y(:,5));
        

        for uu = 1:length(units)

            unit = units(uu);

            if ismember(unit,special_units) && sensitivity(uu,tone_number)
                tone_sensitive = tone_sensitive + 1;
            end
        end
end

%%


load('PNs_05_99.mat')
units_by_rec = PN_units;

tone_sensitive = 0;

for pen = 1:length(units_by_rec)
            
        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/tmp02/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

        stims = unique(type);
        tone_number = find(strcmp(stims,'tone'));
        special_units = units_by_rec{pen,3};
        units = unique(Y(:,3));
        Flist = unique(F0);
        repeats = unique(Y(:,5));
        

        for uu = 1:length(units)

            unit = units(uu);

            if ismember(unit,special_units) && sensitivity(uu,tone_number)
                tone_sensitive = tone_sensitive + 1;
            end
        end
end