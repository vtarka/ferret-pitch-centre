% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, May 2023

% file_names = {'HNs_05_99','TNs_05_99','PNs_05_99'};

figure; 
file_name = 'HNs_05_99';

% units_to_plot = {{'Noah','P08',214},{'Noah','P03',619},{'Ronnie','P05',473}};

units_to_plot = {{'Noah','P08',203},{'Noah','P03',392},{'Ronnie','P05',23}};

% {{'Noah','P02',568},{'Noah','P02',683},{'Ronnie','P13',128}};
% {{'Noah','P08',131},{'Noah','P08',221},{'Dory','P01',370}}
% Updated TNs: {{'Noah','P08',131},{'Dory','P01',370},{'Noah','P08',233}}

% all_units_to_plot = {,...
%     {{'Noah','P08',131},{'Ronnie','P08',231},{'Dory','P01',370}},...
%     {{'Noah','P03',619},{'Noah','P08',214},{'Ronnie','P05',341}}};

stims_to_plot = {'low','high','CT0'}; %{{'low','CT0'},{'high','CT0'},{'low','high','CT0'}};
stims_for_profile = {'low','high','allHarm','CT0','CT5','CT10'};

colors = [0 0 1; 1 0 0; 0 0 0]; %[1 0 0; 0 0 0]; %
ymax = [30 20 30];
% all_colors = {[1 0 0; 0 0 0],[0.3 0 1; 0 0 0],[1 0 0; 0.3 0 1; 0 0 0]};

plot_titles = {'HNs','TNs','PNs'};

mat_struct = load(file_name);
mat_cell = struct2cell(mat_struct);
units_by_rec = mat_cell{1};

window = [0 0.1];


for pen = 1:length(units_to_plot)

    this_unit = units_to_plot{pen};

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' this_unit{1} '/tmp02/Spikes_' this_unit{1} '_' this_unit{2} '_Good_Pitch.mat']);

    subplot(4,2,(pen*2)-1)

    plot_tuning_by_cond(Y,type,F0,this_unit{3},stims_to_plot,zeros(length(stims_to_plot)),this_unit{1},this_unit{2},window,colors);

%     if sp_counter ~= 3
%         xticks([])
%         yticks([])
%         xlabel('')
%         ylabel('')
%     end

%     xticks([])
%     yticks([])
    xlabel('')
    ylabel('')
    sgtitle('')

    ylim([0 ymax(pen)])
    yticks(ymax(pen)/2:ymax(pen)/2:ymax(pen))

    xticks(1:4:17)
    Flist = unique(F0);
    xticklabels(Flist(1:4:17))
    set(gca,'fontsize',18)

    if pen == 1
        xticks(1:4:17)
        Flist = unique(F0);
        xticklabels(Flist(1:4:17))
        xlabel('F0 (Hz)')

        ylh = ylabel('Firing Rate (spks/s)');
%         ylh.Position(1) = ylh.Position(1) - abs(ylh.Position(1) * 0.8);

        set(gca,'fontsize',18)
    end

    rp = get_response_profile(Y,type,F0,this_unit{3},stims_for_profile,window);
    rp = (rp - min(min(rp))) / (max(max(rp)) - min(min(rp)));
    subplot(4,2,pen*2)
    imagesc(rp)
     xticks([])
    yticks([])
    xlabel('')
    ylabel('')
    sgtitle('')

    xticks(1:4:17)
    Flist = unique(F0);
    xticklabels(Flist(1:4:17))
    set(gca,'fontsize',18)

    if pen == 1
        set(gca,'YAxisLocation','right')

        yticks(1:size(rp,1))
        yticklabels({'Low Harm.','High Harm.','Missing F0','Click Train','CT 5% Jitter','CT 10% Jitter'})

        set(gca,'fontsize',18)

        xticks(1:4:17)
        Flist = unique(F0);
        xticklabels(Flist(1:4:17))
        xlabel('F0 (Hz)')
    end

    if pen == 3
        cb = colorbar;
        cb.Position(1) = cb.Position(1) + abs(cb.Position(1) * 0.12);
        cb.Position(4) = cb.Position(4) + abs(cb.Position(4) * 0.8);
        cb.Position(2) = cb.Position(2) + abs(cb.Position(2) * 0.15);

        cb.Ticks = [0 1];
        ylabel(cb,'Norm. Response','FontSize',18,'Rotation',90)
%         cb.Label.Position(1) = cb.Label.Position(1) - abs(cb.Label.Position(1) * 0.2);
    end
   

end


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

% subplot(4,2,7:8)
figure;
C = categorical(loc_list,1:5,{'low A1','high A1','low AAF','high AAF','PPF'});
histogram(C,'FaceColor','k')
ylabel('Number of Units')

ylim([0 20])
set(gca,'fontsize',18)
