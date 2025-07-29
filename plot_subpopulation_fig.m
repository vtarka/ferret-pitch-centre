% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, May 2023


f = figure; 
f.Position = [2623 706 1002 796];

%% HNs
% file_name = '0225HNs';
% null_file = 'HN_nulls';
% units_to_plot = {{'Noah','P02',568},{'Ronnie','P13',128},{'Noah','P02',683}};
% fi = 1;
% stims_to_plot = {'low','high','CT0'};
% stims_for_profile = {'CT0','low','allHarm','CT5'};
% rp_labels = {'click train (CT)','low harm.','missing F0','CT 5% jitter'};
% ex_Ns = [5, 30, 7];
% ymax = [40 80 50];
% colors = [0 0 1; 1 0 0; 0 0 0];
% edges = 0.1:0.1:0.9;

%% TNs
% file_name = '0225TNs';
% null_file = 'TN_nulls';
% units_to_plot = {{'Noah','P03',489},{'Noah','P08',233},{'Dory','P01',370}};
% stims_to_plot = {'high','CT0'};
% stims_for_profile = {'CT0','high','allHarm','CT5'};
% rp_labels = {'click train (CT)','high harm.','missing F0','CT 5% jitter'};
% ex_Ns = [4, 15, 30];
% ymax = [46 80 32];
% edges = 0.1:0.1:0.9;
% fi = 2;
% colors = [1 0 0; 0 0 0];

%% PNs
file_name = '0225PNs';
null_file = 'PN_nulls';
units_to_plot = {{'Noah','P03',619},{'Ronnie','P05',473},{'Ronnie','P13',115}}; 
stims_to_plot = {'low','high','CT0'};
stims_for_profile = {'CT0','low','high','allHarm','CT5'};
fi = 3;
rp_labels = {'click train (CT)','low harm.','high harm.','missing F0','CT 5% jitter'};
ex_Ns = [3, 22, 25];
ymax = [40 96 28];
colors = [0 0 1; 1 0 0; 0 0 0];
edges = 0.1:0.1:0.9;

%%

mat_struct = load(file_name);
mat_cell = struct2cell(mat_struct);
units_by_rec = mat_cell{1};

window = [0 0.1];


for pen = 1:length(units_to_plot)

    this_unit = units_to_plot{pen};

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' this_unit{1} '/Batch0324/Spikes_' this_unit{1} '_' this_unit{2} '_Good_Pitch.mat']);

    subplot(3,3,pen)

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

        ylh = ylabel('firing rate (spks/s)');

        set(gca,'fontsize',18)
        set(gca,'fontname','DejaVu Serif')
    end
    xlabel('F0 (Hz)')

    rp = get_response_profile(Y,type,F0,this_unit{3},stims_for_profile,window);
    rp = (rp - min(min(rp))) / (max(max(rp)) - min(min(rp)));
    subplot(3,3,pen+3)
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
    if pen == 3
        cb = colorbar;
        cb.Position(1) = cb.Position(1) + abs(cb.Position(1) * 0.09);
        cb.Ticks = [0 1];
        ylabel(cb,'Norm. Response','FontSize',10,'Rotation',90)
    end
   

end

%% Plot locations

nUnits = count_units(units_by_rec);
% locs = [1 1 2 1 5 3 5 4 3 3 1 2 3 2 1 1 1 3 4 3];
locs = [1 1 2 1 5 3 5 4 4 1 1 1 1 2 1 1 1 4 2 3];
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
% figure(1);
figure(2);
% subplot(4,2,8)
% subplot(3,3,9)
subplot(1,2,2)
C = categorical(loc_list,1:5,{'low A1','high A1','low AAF','high AAF','PPF'});
histogram(C,'FaceColor','k')
ylabel('Number of Units')

ylim([0 20])
set(gca,'fontsize',18)
set(gca,'fontname','DejaVu Serif')


%% Plot the tuning alignment values

Cs_passing = zeros(nUnits,1);
% stims_for_profile = {'CT0','CT5','CT10','allHarm','low'};
window = [0 0.1];
unit_counter = 1;

% lags = zeros(nUnits,2);
% 
% for pen = 1:length(units_by_rec)
% 
%     load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/tmp02/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);
% 
%     units = units_by_rec{pen,3};
% 
%     for uu = 1:length(units)
%         profile = get_response_profile(Y,type,F0,units(uu),stims_for_profile,window);
% 
%         profile = profile(:,3:17); % eliminate the first two frequencies because they weren't presented for low
% 
%         real_corr = get_avg_pairwise_corr(profile);
%         Cs(unit_counter) = real_corr;
% 
%         [r,l] = xcorr(profile(2,:),profile(4,:));
%         [max_r,max_r_idx] = max(r);
%         lags(unit_counter,1) = max_r;
%         lags(unit_counter,2) = l(max_r_idx);
% 
%         unit_counter = unit_counter +1;
%     end
% end
% 
% figure;
% histogram(lags(:,2))
% set(gca,'fontsize',20)
% xlim([-10 10])
% xticks([-8 -4 -2 0 2 4 8])
% xticklabels([-2 -1 -0.5 0 0.5 1 2])
% % ylim([0 9])
% yticks(0:3:12)
% 
% xlabel('Best Lag (octaves)')
% ylabel('Neuron Count')
% 
% figure(1);
% 
% % subplot(4,2,7)
% subplot(3,3,7)
% histogram(Cs,edges,'facecolor','k')
% set(gca,'fontsize',18)
% set(gca,'fontname','DejaVu Serif')
% xlim([0.1 0.9])
% 
% if fi == 1
%     yticks(0:4:12)
%     ylim([0 13])
% elseif fi == 2
%     yticks(0:4:12)
%     ylim([0 10])
% else
%     yticks(0:5:15)
%     ylim([0 15])
% end
% 
% xlabel('Average Correlation')
% ylabel('Neurons')
% xticks(0.1:0.2:0.9)


nDiscarded = count_units(units_by_rec,4);
Cs_passing = zeros(nUnits,1);
Cs_discarded = zeros(nDiscarded,1);
% stims_for_profile = {'CT0','CT5','CT10','allHarm','low'};
window = [0 0.1];
unit_counter = 1;
for pen = 1:length(units_by_rec)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/Batch0324/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

    units = units_by_rec{pen,3};

    for uu = 1:length(units)
        profile = get_response_profile(Y,type,F0,units(uu),stims_for_profile,window);

        profile = profile(:,3:17); % eliminate the first two frequencies because they weren't presented for low

        real_corr = get_avg_pairwise_corr(profile);
        Cs_passing(unit_counter) = real_corr;

%         if unit_counter == ex_Ns(1)
%             figure(1);
% %             subplot(4,2,7)
%             subplot(3,3,7)
%             hold on;
% %             scatter(real_corr,0,100,[0,1,0],'filled','markerfacealpha',1,'markerfacecolor',[0,1,0],'markeredgecolor','k')
%             scatter(real_corr,0,75,'filled','k','LineWidth',2,'markerfacealpha',0.4,'markeredgecolor','k')
%         elseif unit_counter == ex_Ns(2)
%             figure(1);
% %             subplot(4,2,7)
%             subplot(3,3,7)
%             hold on;
% %             scatter(real_corr,0,100,[1,0,1],'filled','markerfacealpha',1,'markerfacecolor',[1,0,1],'markeredgecolor','k')
%             scatter(real_corr,0,75,'filled','k','d','LineWidth',2,'markerfacealpha',0.4,'markeredgecolor','k')
%         elseif unit_counter == ex_Ns(3)
%             figure(1);
% %             subplot(4,2,7)
%             subplot(3,3,7)
%             hold on;
% %             scatter(real_corr,0,100,[1,1,0],'filled','markerfacealpha',1,'markerfacecolor',[1,1,0],'markeredgecolor','k')
%             scatter(real_corr,0,250,'filled','k','pentagram','LineWidth',2,'markerfacealpha',0.4,'markeredgecolor','k')
%         end
        unit_counter = unit_counter +1;
    end

    discarded_units = units_by_rec{pen,4};
    for uu = 1:length(discarded_units)
        profile = get_response_profile(Y,type,F0,discarded_units(uu),stims_for_profile,window);

        profile = profile(:,3:17); % eliminate the first two frequencies because they weren't presented for low

        real_corr = get_avg_pairwise_corr(profile);
        Cs_discarded(uu) = real_corr;
       
    end
end

hold off;


set(findobj(gcf,'type','axes'),'FontName','DejaVu Serif','FontSize',10);

figure(2);
subplot(1,2,1)
hold on
histogram(Cs_passing,'binwidth',0.05,'facecolor','black','facealpha',0.5,'linewidth',2)
histogram(Cs_discarded,'binwidth',0.05,'facecolor','white','facealpha',0.5,'linewidth',2)

legend({'harmonicity neurons','neurons with non-significant tuning correlation'},'fontsize',14)
xlabel('average tuning correlation')
ylabel('# of neurons')

set(gca,'fontsize',20)


%% plot null distributions for example neurons

load(null_file)
figure(1);

% units_to_plot = {{'Noah','P02',568},{'Ronnie','P13',128},{'Noah','P02',683}};

for i = 1:3

    unit_to_plot = units_to_plot{i};

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit_to_plot{1} '/Batch0324/Spikes_' unit_to_plot{1} '_' unit_to_plot{2} '_Good_Pitch.mat']);

    this_unit = unit_to_plot{3};

    rp = get_response_profile(Y,type,F0,this_unit,stims_for_profile,window);
    avg_corr = get_avg_pairwise_corr(rp);
    
    subplot(3,3,i+6)
    hold on

    xline(prctile(nullDistributions(ex_Ns(i),:),95),'linewidth',5,'linestyle',':','color','k')
    xline(avg_corr,'linewidth',5,'color','k')
    histogram(nullDistributions(ex_Ns(i),:),'FaceColor','black','facealpha',0.25)

    xlim([-0.3 0.8])

    

%     set(gca,'fontsize',20)

    if i == 1
        xlabel('average correlation')
        ylabel('shuffled simulations')
    end

    if i == 3
        legend({'95th %ile','observed corr.','null distribution'},'fontsize',10)
    end

    xlim([-0.3 0.8])
    xticks(-0.2:0.2:0.9)
    yticks(0:200:800)

    set(gca,'fontsize',18)
    set(gca,'fontname','DejaVu Serif')

end

