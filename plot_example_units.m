%% Plot 3 examples of units in each population
% DEPENDENCIES:
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, April 2023


figure; 
sps = [1 4 7; 2 5 8; 3 6 9]';
sp_counter = 1;

file_names = {'HNs_05_99','TNs_05_99','PNs_05_99'};

all_units_to_plot = {{{'Noah','P02',568},{'Noah','P02',683},{'Ronnie','P13',128}},...
    {{'Noah','P08',131},{'Ronnie','P08',231},{'Dory','P01',370}},...
    {{'Noah','P03',619},{'Noah','P08',214},{'Ronnie','P05',341}}};

stims_to_plot = {{'low','CT0'},{'high','CT0'},{'low','high','CT0'}};

all_colors = {[1 0 0; 0 0 0],[0.3 0 1; 0 0 0],[1 0 0; 0.3 0 1; 0 0 0]};

plot_titles = {'HNs','TNs','PNs'};

window = [0 0.1];

for ff = 1:length(file_names)

    mat_struct = load(file_names{ff});
    mat_cell = struct2cell(mat_struct);
    units_by_rec = mat_cell{1};

    colors = all_colors{ff};

    units_to_plot = all_units_to_plot{ff};
    
    for pen = 1:length(units_to_plot)

        this_unit = units_to_plot{pen};

        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' this_unit{1} '/tmp02/Spikes_' this_unit{1} '_' this_unit{2} '_Good_Pitch.mat']);

        subplot(3,3,sps(sp_counter))

        plot_tuning_by_cond(Y,type,F0,this_unit{3},stims_to_plot{ff},zeros(length(stims_to_plot{ff})),this_unit{1},this_unit{2},window,colors)

        if sp_counter ~= 3
            xticks([])
            yticks([])
            xlabel('')
            ylabel('')
        end

        sgtitle('')

        sp_counter = sp_counter + 1;

    end
end
