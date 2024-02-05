% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, May 2023

% file_names = {'HNs_05_99','TNs_05_99','PNs_05_99'};

f = figure; 
f.Position = [1654 433 1900 1333];

% units_to_plot = {{'Noah','P08',214},{'Noah','P03',619},{'Ronnie','P05',473}};

%% HNs
% file_name = 'HNs_05_99_3';
% units_to_plot = {{'Noah','P02',568},{'Ronnie','P13',128},{'Noah','P02',683}};
% fi = 1;
% stims_to_plot = {'low','high','CT0'};
% stims_for_profile = {'low','allHarm','CT0','CT5','CT10'};
% rp_labels = {'Low Harm.','Missing F0','Click Train','CT 5% Jitter','CT 10% Jitter'};
% ex_Ns = [5, 17, 7];
% ymax = [40 80 50];
% colors = [0 0 1; 1 0 0; 0 0 0];
% edges = 0.2:0.1:0.7;

%% TNs
% file_name = 'TN_units_05_99_2';
% units_to_plot = {{'Noah','P03',489},{'Noah','P08',233},{'Dory','P01',370}};
% stims_to_plot = {'high','CT0'};
% stims_for_profile = {'high','allHarm','CT0','CT5','CT10'};
% rp_labels = {'High Harm.','Missing F0','Click Train','CT 5% Jitter','CT 10% Jitter'};
% ex_Ns = [9, 21, 43];
% ymax = [46 80 32];
% edges = 0.2:0.1:0.9;
% fi = 2;
% colors = [1 0 0; 0 0 0];

%% PNs
file_name = 'PNs_05_99_3';
units_to_plot = {{'Noah','P03',619},{'Ronnie','P05',473},{'Ronnie','P13',115}}; 
stims_to_plot = {'low','high','CT0'};
stims_for_profile = {'low','high','allHarm','CT0','CT5','CT10'};
fi = 3;
rp_labels = {'Low Harm.','High Harm.','Missing F0','Click Train','CT 5% Jitter','CT 10% Jitter'};
ex_Ns = [3, 22, 25];
ymax = [40 82 28];
colors = [0 0 1; 1 0 0; 0 0 0];
edges = 0.1:0.1:0.7;
subtype = 'PN';

%% Plot 3 Example Neurons

mat_struct = load(file_name);
mat_cell = struct2cell(mat_struct);
units_by_rec = mat_cell{1};

window = [0 0.1];


for pen = 1:length(units_to_plot)

    this_unit = units_to_plot{pen};

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' this_unit{1} '/tmp02/Spikes_' this_unit{1} '_' this_unit{2} '_Good_Pitch.mat']);

    figure;

    plot_tuning_by_cond(Y,type,F0,this_unit{3},stims_to_plot,zeros(length(stims_to_plot)),this_unit{1},this_unit{2},window,colors);

    xlabel('')
    ylabel('')
    sgtitle('')
    set(gca,'fontname','DejaVu Serif')

    ylim([0 ymax(pen)])
    yticks(ymax(pen)/2:ymax(pen)/2:ymax(pen))

    xticks(1:4:17)
    Flist = unique(F0);
    xticklabels(Flist(1:4:17))
    set(gca,'fontsize',18)
    set(gca,'fontname','DejaVu Serif')

    if pen == 1
        xticks(1:4:17)
        Flist = unique(F0);
        xticklabels(Flist(1:4:17))
        xlabel('F0 (Hz)')

        ylh = ylabel('Firing Rate (spks/s)');

        set(gca,'fontsize',18)
        set(gca,'fontname','DejaVu Serif')
    end

    xlabel('F0 (Hz)')

    fname = sprintf('%s_%s_%d',subtype,'tuning',pen);
    saveas(gcf,fname,'fig')

    rp = get_response_profile(Y,type,F0,this_unit{3},stims_for_profile,window);
    rp = (rp - min(min(rp))) / (max(max(rp)) - min(min(rp)));

    figure;
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
    set(gca,'fontname','DejaVu Serif')

    if pen == 1
%         set(gca,'YAxisLocation','right')

        yticks(1:size(rp,1))
        yticklabels(rp_labels)

        set(gca,'fontsize',18)
        set(gca,'fontname','DejaVu Serif')

        xticks(1:4:17)
        Flist = unique(F0);
        xticklabels(Flist(1:4:17))
        xlabel('F0 (Hz)')
    end

    xlabel('F0 (Hz)')

    fname = sprintf('%s_%s_%d',subtype,'rp',pen);
    saveas(gcf,fname,'fig')

    if pen == 3

        figure;
        cb = colorbar;
%         cb.Position(1) = cb.Position(1) + abs(cb.Position(1) * 0.08);
%         cb.Position(4) = cb.Position(4) + abs(cb.Position(4) * 0.8);
%         cb.Position(2) = cb.Position(2) + abs(cb.Position(2) * 0.15);

        cb.Ticks = [0 1];
        ylabel(cb,'Norm. Response','FontSize',18,'Rotation',90)
%         cb.Label.Position(1) = cb.Label.Position(1) - abs(cb.Label.Position(1) * 0.2);

        saveas(gcf,'colorbar.fig')
    end
   
    

end

%% Plot locations

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
% subplot(4,2,8)
% subplot(3,3,9)
C = categorical(loc_list,1:5,{'low A1','high A1','low AAF','high AAF','PPF'});
histogram(C,'FaceColor','k')
ylabel('Number of Units')

ylim([0 20])
set(gca,'fontsize',18)
set(gca,'fontname','DejaVu Serif')

fname = sprintf('%s_locations',subtype);
saveas(gcf,fname,'fig')


%% Plot the tuning alignment values

Cs = zeros(nUnits,1);
window = [0 0.1];
unit_counter = 1;

for pen = 1:length(units_by_rec)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/final/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

    units = units_by_rec{pen,3};

    for uu = 1:length(units)
        profile = get_response_profile(Y,type,F0,units(uu),stims_for_profile,window);

        profile = profile(:,3:17); % eliminate the first two frequencies because they weren't presented for low

        real_corr = get_avg_pairwise_corr(profile);
        Cs(unit_counter) = real_corr;
        unit_counter = unit_counter +1;
    end
end

figure;
histogram(Cs,edges,'facecolor','k')
set(gca,'fontsize',18)
set(gca,'fontname','DejaVu Serif')

if fi == 1
    yticks(0:4:12)
    ylim([0 12])
elseif fi == 2
    yticks(0:4:12)
    ylim([0 12])
else
    yticks(0:5:15)
    ylim([0 15])
end

xlabel('Average Correlation')
ylabel('Neurons')
xticks(0.1:0.2:0.9)


Cs = zeros(nUnits,1);
% stims_for_profile = {'CT0','CT5','CT10','allHarm','low'};
window = [0 0.1];
unit_counter = 1;
for pen = 1:length(units_by_rec)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/final/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

    units = units_by_rec{pen,3};

    for uu = 1:length(units)
        profile = get_response_profile(Y,type,F0,units(uu),stims_for_profile,window);

        profile = profile(:,3:17); % eliminate the first two frequencies because they weren't presented for low

        real_corr = get_avg_pairwise_corr(profile);
        Cs(unit_counter) = real_corr;

        if unit_counter == ex_Ns(1)
            hold on;
            scatter(real_corr,0,150,'filled','k','LineWidth',2,'markerfacealpha',0.4,'markeredgecolor','k')
        elseif unit_counter == ex_Ns(2)
            hold on;
            scatter(real_corr,0,150,'filled','k','d','LineWidth',2,'markerfacealpha',0.4,'markeredgecolor','k')
        elseif unit_counter == ex_Ns(3)
            hold on;
            scatter(real_corr,0,350,'filled','k','pentagram','LineWidth',2,'markerfacealpha',0.4,'markeredgecolor','k')
        end
        unit_counter = unit_counter +1;
    end
end

hold off;


set(findobj(gcf,'type','axes'),'FontName','DejaVu Serif','FontSize',18);

fname = sprintf('%s_tuning_alignment',subtype);

saveas(gcf,fname,'fig')

