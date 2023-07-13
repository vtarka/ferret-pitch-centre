%% Plot the locations of the provided population (temporal, harmonicity, or pitch)
% DEPENDENCIES:
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, April 2023

figure; 

file_names = {'HNs_05_99','TNs_05_99','PNs_05_99'};
plot_titles = {'HNs','TNs','PNs'};

for ff = 1:length(file_names)

    mat_struct = load(file_names{ff});
    mat_cell = struct2cell(mat_struct);
    units_by_rec = mat_cell{1};
    
    nUnits = count_units(units_by_rec);
    locs = [1 1 2 1 5 3 5 4 3 3 1 2 3 2 1 1 1 3 4 3];
    loc_list = zeros(nUnits,1);
    unit_counter = 1;
    
    for pen = 1:length(units_by_rec)
    
        units = units_by_rec{pen,3};
    
        for uu = 1:length(units)
            loc_list(unit_counter) = locs(pen);
            unit_counter = unit_counter+1;
        end
    end
    
    subplot(1,3,ff)
    C = categorical(loc_list,1:4,{'low A1','high A1','low AAF','high AAF'});
    histogram(C,'FaceColor','k')
    ylabel('Number of Units')
    set(gca,'fontsize',22)

    ylim([0 20])

    title(plot_titles{ff})
end
