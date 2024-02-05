
figure; 

file_name = 'TN_units_f2';
units_to_plot = {{'Noah','P03',634},{'Noah','P08',237},{'Noah','P08',221},{'Ronnie','P05',306}}; 
stims_to_plot = {'high','alt','rand'};
% stims_for_profile = {'low','high','allHarm','CT0','CT5','CT10'};
% fi = 3;
% rp_labels = {'Low Harm.','High Harm.','Missing F0','Click Train','CT 5% Jitter','CT 10% Jitter'};
% ex_Ns = [3, 22, 25];
% ymax = [40 82 28];
colors = [0 0 0; 0.85 0.325 0.098;0.466 0.674 0.188];
edges = 0.1:0.1:0.7;
ymax = [36 60 100 42];

mat_struct = load(file_name);
mat_cell = struct2cell(mat_struct);
units_by_rec = mat_cell{1};

window = [0 0.1];

SPs = [3 5 4 6];


for pen = 1:length(units_to_plot)

    this_unit = units_to_plot{pen};

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' this_unit{1} '/tmp02/Spikes_' this_unit{1} '_' this_unit{2} '_Good_Pitch.mat']);

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

    xlabel('F0 (Hz)')

    if pen == 1
        ylabel('Firing Rate (spikes/s)')
    end
end


subplot(3,2,1)
bar([7 38-7],'k')
set(gca,'fontsize',18)
xticks(1:2)
xticklabels({'Sensitive','Not Sensitive'})
yticks(0:10:30)
ylim([0 32])
ylabel('Temporal Neurons')