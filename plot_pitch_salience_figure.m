
figure; 

% file_name = 'TNs_05_99_3';
units_to_plot = {{'Ronnie','P05',98},{'Noah','P02',683},{'Ronnie','P04',219},{'Ronnie','P05',410}}; 
stims_to_plot = {'CT0','CT5','CT10','CT20','CT40'};
colors = [1 0 0; 1 0 1; 0 0 0; 0 1 1; 0 0 1];
edges = 0.1:0.1:0.7;
ymax = [90 40 20 44];

% mat_struct = load(file_name);
% mat_cell = struct2cell(mat_struct);
% units_by_rec = mat_cell{1};
window = [0 0.1];

SPs = [3 5 4 6];


for pen = 1:length(units_to_plot)

    this_unit = units_to_plot{pen};

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' this_unit{1} '/final/Spikes_' this_unit{1} '_' this_unit{2} '_Good_Pitch.mat']);

%     subplot(4,2,(pen*2)-1)
    subplot(3,2, SPs(pen))

    plot_tuning_by_cond(Y,type,F0,this_unit{3},stims_to_plot,zeros(length(stims_to_plot)),this_unit{1},this_unit{2},window,colors);

    xlabel('')
    ylabel('')
    sgtitle('')

    ylim([0 ymax(pen)])
    yticks(ymax(pen)/2:ymax(pen)/2:ymax(pen))

    xticks(1:4:17)
    Flist = unique(F0);
    xticklabels(Flist(1:4:17))
    set(gca,'fontsize',18)
    set(gca,'fontname','DejaVu Serif')

    xlabel('F0 (Hz)')

    if pen == 1 || pen == 2
        ylabel('Firing Rate (spikes/s)')
    end
end