%% Plot all the response profiles in a single population
% DEPENDENCIES: get_response_profile.m
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, April 2023


mat_struct = load('TNs_01_99.mat');
mat_cell = struct2cell(mat_struct);
units_by_rec = mat_cell{1};

stims = {'CT0','CT5','CT10','allHarm','high'};

window = [0 0.1]; % in seconds

figure;

for pen = 1:length(units_by_rec)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/tmp02/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

    units = units_by_rec{pen,3};

    for uu = 1:length(units)
        unit = units(uu);

        nexttile;
        imagesc(get_response_profile(Y,type,F0,unit,stims,window))
        xticks([])
        yticks([])

    end
end