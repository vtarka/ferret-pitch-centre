%% Analyse depths of neurons for Noah recordings

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08'};

all_depths = cell(8,1);


% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/final/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    depth_file = ['/media/veronica/Kat Data/Noah/For_analysis/' Pens{ap} '/' Pens{ap} '-pitch70dB2020_g0/cluster_info.tsv'];

    t = readtable(depth_file,'FileType','text','Delimiter','\t');

    depths = [t.cluster_id t.depth];

    units = unique(Y(:,3));

    units_logical = ismember(depths(:,1),units);

    unit_depths = depths(units_logical,:);

    all_depths{ap} = unit_depths;

end


%% plot depth distributions


load('TN_units_new_05.mat')

load('HN_units_new_05.mat')


figure;

titles = {'Noah 1 (low A1)','Noah 2 (low A1)','Noah 3 (high A1)','Noah 4 (low A1)','Noah 5 (low PPF)','Noah 6 (low AAF)','Noah 7 (low PPF)','Noah 8 (high AAF)'};

for i = 1:8

    subplot(2,4,i)
    hold on

    depths = all_depths{i};

    histogram(depths(:,2),0:500:4000,'Normalization','probability','FaceColor','white','FaceAlpha',0.2,'linewidth',2)

    
    % HNs

    HNs = HN_units{i,3};

    HN_depths_logical = ismember(depths(:,1),HNs);

    HN_depths = depths(HN_depths_logical,2);

    histogram(HN_depths,0:500:4000,'Normalization','probability','FaceColor','blue','facealpha',0.3)


    % TNs
    
    TNs = TN_units{i,3};

    TN_depths_logical = ismember(depths(:,1),TNs);

    TN_depths = depths(TN_depths_logical,2);

    histogram(TN_depths,0:500:4000,'Normalization','probability','FaceColor','red','facealpha',0.3)
    
    
    set(gca,'fontsize',20)

    title(titles{i})

end

subplot(2,4,5)

legend({'All Neurons','Harmonicity Neurons','Temporal Neurons'},'fontsize',14)

subplot(2,4,1)
xlabel('Depth')
ylabel('% of Neurons')

%% plot HN depth distributions

load('HN_units_new_05.mat')
figure;

for ap = 1:8

    HNs = HN_units{ap,3};

    depths = all_depths{ap};

    HN_depths_logical = ismember(depths(:,1),TNs);

    HN_depths = depths(HN_depths_logical,2);

    subplot(2,4,ap)

    histogram(HN_depths,0:500:4000)
    set(gca,'fontsize',20)
    title(titles{ap})

end

sgtitle('Harmonicity Neurons')
subplot(2,4,1)
xlabel('Depth')
ylabel('# of Neurons')

%% plot TN depth distributions

load('TN_units_new_05.mat')
figure;

for ap = 1:8

    TNs = TN_units{ap,3};

    depths = all_depths{ap};

    TN_depths_logical = ismember(depths(:,1),TNs);

    TN_depths = depths(TN_depths_logical,2);

    subplot(2,4,ap)

    histogram(TN_depths,0:500:4000)
    set(gca,'fontsize',20)
    title(titles{ap})

end

sgtitle('Temporal Neurons')

subplot(2,4,1)
xlabel('Depth')
ylabel('# of Neurons')