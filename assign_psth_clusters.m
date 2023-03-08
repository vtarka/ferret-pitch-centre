
%% Clusters neurons based on recording-averaged PSTHs
% OUTPUT: saves an additional variable to the spiketime files called
% 'clusters' nUnits x 1 with the cluster assignment for each unit
% DEPENDENCIES:
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, March 2023

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';


% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

windows = [0:10:630; 10:10:640]'/1000; % 10 ms windows to sum spiketimes over to build the PSTH

% initialize the space to save the data we'll put into k-means
% every unit has one row, which is an average of all of its responses
nUnits = 1311;
all_unit_psth = zeros(nUnits,length(windows)); 
pen_labels = zeros(nUnits,1);
loc_labels = zeros(nUnits,1);

% functional locations of each penetration (as listed in variable 'Pens')
%   low A1      high A1       low AAF      high AAF      PPF
%     1            2             3           4            5
locs = [1 1 2 1 5 3 5 4 3 3 1 2 3 2 1 1 1 3 4 3]; 

uCounter = 1;

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp01/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    units = unique(Y(:,3));

    for uu = 1:length(units)

        unitSpikes = Y(Y(:,3)==units(uu),:);

        trials = unique(unitSpikes(:,6));
        nTrials = length(trials);

        nSpikes = zeros(nTrials,length(windows));

        for ww = 1:length(windows)

            win_start = windows(ww,1);
            win_end = windows(ww,2);

            [r,c] = find(unitSpikes(:,2) > win_start & unitSpikes(:,2) <= win_end);
            trials_with_spikes = unique(Y(r,6));

            for tt = 1:nTrials

                if ismember(trials(tt),trials_with_spikes)

                    spikes = unitSpikes(:,2) > win_start & unitSpikes(:,2) <= win_end & unitSpikes(:,6) == trials(tt);
                    nSpikes(tt,ww) = sum(spikes);

                end

            end % ends window loop

        end % ends trial loop

        all_unit_psth(uCounter,:) = mean(nSpikes);
        pen_labels(uCounter) = ap;
        loc_labels(uCounter) = locs(ap);

        uCounter = uCounter + 1;

    end % ends unit loop
end % ends recording loop

%% Adjust Noah's times

NoahUnits = all_unit_psth(pen_labels<9,:);

for nu = 1:length(NoahUnits)

    offset200 = NoahUnits(nu,30:59);
    NoahUnits(20:49) = offset200;

end

all_unit_psth(pen_labels<9,:) = NoahUnits;
all_unit_psth(:,50:end) = [];

%% Eliminate PSTHs that are all zeros

active_units = zeros(size(all_unit_psth,1),1);

% for every unit
for uu = 1:length(all_unit_psth)

    % check if there were no spikes whatsoever for this unit
    if ~isempty(find(all_unit_psth(uu,:),1))
        active_units(uu) = 1;
    end

end

active_all_unit_psth = all_unit_psth(active_units==1,:);
active_pen_labels = pen_labels(active_units==1);
active_loc_labels = loc_labels(active_units==1);

%% Normalize the PSTHs

norm_active_psths = zeros(size(active_all_unit_psth));

for uu = 1:length(active_all_unit_psth)

    mean_subtracted = active_all_unit_psth(uu,:) - mean(active_all_unit_psth(uu,:));
    std_normalized = mean_subtracted / std(active_all_unit_psth(uu,:));
    norm_active_psths(uu,:) = std_normalized;

end

within_cluster_var = zeros(25,1);
for i = 1:25

    k = i;
    [~,~,sumd,~] = kmeans(norm_active_psths,k);
    within_cluster_var(i) = mean(sumd);

end

%% Generate the clusters based on the optimal determined from within_cluster_var

k = 5;
[idx,C,sumd,D] = kmeans(norm_active_psths,k);

%% Plot the clusters

figure;

cluster_sizes = zeros(k,1);

loc_frequencies = zeros(5,k);
pen_frequencies = zeros(20,k);

for i = 1:k
    
    thisCluster = norm_active_psths(idx==i,:);
    cluster_sizes(i) = size(thisCluster,1);

    cluster_locs = active_loc_labels(idx==i);

    [nLoc,edgesLoc] = histcounts(cluster_locs,5);

    norm_locs = zeros(size(nLoc));
    for ll = 1:length(edgesLoc)-1
        norm_locs(ll) = nLoc(ll)/length(find(active_loc_labels==ll));
    end
    loc_frequencies(:,i) = norm_locs;

    cluster_pens = active_pen_labels(idx==i);

    [nPen,edgesPen] = histcounts(cluster_pens,20);

    norm_pens = zeros(size(nPen));
    for ll = 1:length(edgesPen)-1
        norm_pens(ll) = nPen(ll)/length(find(active_pen_labels==ll));
    end
    pen_frequencies(:,i) = norm_pens;

    if ~isempty(thisCluster)

        subplot(3,2,i); errorbar(mean(thisCluster),std(thisCluster),'LineWidth',3)
        title(sprintf('Cluster %d',i))
        xticks(0:15:60)
        xticklabels(0:150:600)
        axis tight

        set(gca,'fontsize',22)
    end

end

ylabel('Normalized Firing Rate')
xlabel('ms after stimulus onset')

figure;

bar(cluster_sizes)
title('Number of Neurons per Cluster')
yticks(0:100:400)
set(gca,'Fontsize',22)

%% Save the cluster assignments to the spiketime variables

uCounter = 1;
active_uCounter = 1;

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    units = unique(Y(:,3));
    clusters = zeros(length(units),1); % initialize space to save the cluster information

    for uu = 1:length(units)

        if active_units(uCounter) == 1
            clusters(uu) = idx(active_uCounter);
            active_uCounter = active_uCounter + 1;
        end

        uCounter = uCounter + 1;
    end

    save(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp01/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat'],...
       'Y','type','F0','sensitivity','clusters')
end
