%% Script to find clusters of neurons by PSTH shape
% DEPENDENCES:
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, February 2023

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';


% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

windows = [0:10:590; 10:10:600]'/1000;

% initialize the space to save the data we'll put into k-means
% every unit has one row, which is an average of all of its responses
all_unit_psth = zeros(1311,length(windows)); 
uCounter = 1;

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp01/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    units = unique(Y(:,3));

    for uu = 1:length(units)

        if ~isempty(find(active(uu,:,1), 1))

            unitSpikes = Y(Y(:,3)==units(uu),:);

            trials = unique(unitSpikes(:,6));
            nTrials = length(trials);

            nSpikes = zeros(nTrials,length(windows));

            for tt = 1:nTrials

                for ww = 1:length(windows)
            
                win_start = windows(ww,1);
                win_end = windows(ww,2);

                spikes = unitSpikes(:,2) >= win_start & unitSpikes(:,2) < win_end & unitSpikes(:,6) == trials(tt);
                nSpikes(tt,ww) = sum(spikes);

                end % ends window loop

            end % ends trial loop

            all_unit_psth(uCounter,:) = mean(nSpikes);
            
        end % ends active if condition

        uCounter = uCounter + 1;

    end % ends unit loop
end % ends recording loop


%% Label which penetrations the units came from

uCounter = 1;
penLabels = zeros(length(all_unit_psth),1);

for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp01/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    units = unique(Y(:,3));

    for uu = 1:length(units)
        penLabels(uCounter) = ap;
        uCounter = uCounter + 1;
    end
end

%% Remove extra stimulus time from Noah units

NoahUnits = all_unit_psth(penLabels<9,:);

for nu = 1:length(NoahUnits)

    postStim100 = NoahUnits(nu,30:39);
    replacement100 = NoahUnits(nu,40:49);
    NoahUnits(nu,20:29) = postStim100;
    NoahUnits(nu,30:39) = replacement100;
end

all_unit_psth(penLabels<9,:) = NoahUnits;

%% Normalize and cluster the resulting matrix

active_all_unit_psth = zeros(size(all_unit_psth));
activePenLabels = zeros(length(penLabels),1);
emptyCounter = 0;
activeCounter = 1;

for uu = 1:length(all_unit_psth)

    if ~isempty(find(all_unit_psth(uu,:),1))
        active_all_unit_psth(activeCounter,:) = all_unit_psth(uu,:);
        activePenLabels(activeCounter) = penLabels(uu);
        activeCounter = activeCounter + 1;
    else
        emptyCounter = emptyCounter + 1;
    end

end

active_all_unit_psth(end-emptyCounter-1:end,:) = [];
activePenLabels(end-emptyCounter-1:end) = [];

norm_active_psths = zeros(size(active_all_unit_psth));
for uu = 1:length(active_all_unit_psth)
        
    mean_subtracted = active_all_unit_psth(uu,:) - mean(active_all_unit_psth(uu,:));
    range_normalized = mean_subtracted / range(active_all_unit_psth(uu,:));
    norm_active_psths(uu,:) = range_normalized;

end

within_cluster_var = zeros(25,1);
for i = 1:25

    k = i;
    [~,~,sumd,~] = kmeans(norm_active_psths,k);
    within_cluster_var(i) = mean(sumd);

end

%%

k = 5;
[idx,C,sumd,D] = kmeans(norm_active_psths,k);

figure;

clusterSizes = zeros(k,1);
for i = 1:k
    
    thisCluster = norm_active_psths(idx==i,:);
    clusterSizes(i) = size(thisCluster,1);

    if ~isempty(thisCluster)

        subplot(3,2,i); errorbar(mean(thisCluster),std(thisCluster),'LineWidth',3)
        title(sprintf('Cluster %d',i))
        xticks(0:15:60)
        xticklabels(0:150:600)
        yticks(-0.5:0.3:1.5)

        set(gca,'fontsize',22)
    end

end

ylabel('Normalized Firing Rate')
xlabel('ms after stimulus onset')
subplot(3,2,6); bar(clusterSizes)
title('Number of Neurons per Cluster')
yticks(0:100:400)
set(gca,'Fontsize',22)

%%
figure;

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

incl = [4 9 10 12 14];
sp_counter = 1;
for ap = 1:length(Animals)

    if ismember(ap,incl)
        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp01/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);
    
        subplot(2,3,sp_counter) 
        histogram(Y(Y(:,2)~=0,2))
        xlabel('Spike time (s after stimulus)')
        sp_counter = sp_counter + 1;
    end
end
%     units = unique(Y(:,3));
% 
%     for uu = 1:length(units)
%         unitSpikes = Y(Y(:,3)==units(uu),:);
%         clf
%         histogram(unitSpikes(:,2),0:.01:max(unitSpikes(:,2)))
%         pause
%     end
% end
