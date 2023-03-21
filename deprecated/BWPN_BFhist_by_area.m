
%% Looking at PN best frequency histogram per area
% Veronica Tarka
% veronica.tarka@dpag.ox.ac.uk
% January 2023

%%
% Low / Mid A1 (Noah 2, Ronnie 8, Ronnie 13)

PN_CFs = [];
nPN_CFs = [];

recording_idx = find(strcmp(PN_units(:,1),'Noah') & strcmp(PN_units(:,2),'P02'));
pen = recording_idx;
load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' PN_units{recording_idx,1}...
    '/tmp/Spikes_' PN_units{recording_idx,1} '_' PN_units{recording_idx,2} '_Good_Pitch.mat']);

stims = unique(type);
Flist = unique(F0);
repeats = unique(Y(:,5));
allUnits = unique(Y(:,3));
PNUnits = PN_units{pen,3};

[~,PNUnit_IDXs] = ismember(PNUnits,allUnits);
PN_tone_CFs = BFs(PNUnit_IDXs,13);
PN_CFs = [PN_CFs; PN_tone_CFs];
nPN_tone_CFs = BFs(:,13);
nPN_tone_CFs(PNUnit_IDXs) = [];
nPN_tone_CFs(nPN_tone_CFs==0) = [];
nPN_CFs = [nPN_CFs; nPN_tone_CFs];

recording_idx = find(strcmp(PN_units(:,1),'Ronnie') & strcmp(PN_units(:,2),'P08'));
pen = recording_idx;
load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' PN_units{recording_idx,1}...
    '/tmp/Spikes_' PN_units{recording_idx,1} '_' PN_units{recording_idx,2} '_Good_Pitch.mat']);

stims = unique(type);
Flist = unique(F0);
repeats = unique(Y(:,5));
allUnits = unique(Y(:,3));
PNUnits = PN_units{pen,3};

[~,PNUnit_IDXs] = ismember(PNUnits,allUnits);
PN_tone_CFs = BFs(PNUnit_IDXs,13);
PN_CFs = [PN_CFs; PN_tone_CFs];
nPN_tone_CFs = BFs(:,13);
nPN_tone_CFs(PNUnit_IDXs) = [];
nPN_tone_CFs(nPN_tone_CFs==0) = [];
nPN_CFs = [nPN_CFs; nPN_tone_CFs];

% recording_idx = find(strcmp(PN_units(:,1),'Ronnie') & strcmp(PN_units(:,2),'P13'));
% pen = recording_idx;
% load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' PN_units{recording_idx,1}...
%     '/tmp/Spikes_' PN_units{recording_idx,1} '_' PN_units{recording_idx,2} '_Good_Pitch.mat']);
% 
% stims = unique(type);
% Flist = unique(F0);
% repeats = unique(Y(:,5));
% allUnits = unique(Y(:,3));
% PNUnits = PN_units{pen,3};
% 
% [~,PNUnit_IDXs] = ismember(PNUnits,allUnits);
% PN_tone_CFs = BFs(PNUnit_IDXs,13);
% PN_CFs = [PN_CFs; PN_tone_CFs];
% nPN_tone_CFs = BFs(:,13);
% nPN_tone_CFs(PNUnit_IDXs) = [];
% nPN_tone_CFs(nPN_tone_CFs==0) = [];
% nPN_CFs = [nPN_CFs; nPN_tone_CFs];


figure;
[PN_N,PN_edges] = histcounts(PN_CFs,1:18);
PN_CFs_norm = PN_N/length(PN_CFs);

[nPN_N,~] = histcounts(nPN_CFs,1:18);
nPN_CFs_norm = nPN_N/length(nPN_CFs);
colors = colormap(jet(2));
histogram('Categories',string(Flist),'BinCounts',PN_CFs_norm,'FaceColor',colors(1,:),'FaceAlpha',0.5)
hold on
histogram('Categories',string(Flist),'BinCounts',nPN_CFs_norm,'FaceColor',colors(2,:),'FaceAlpha',0.5)
ylim([0 0.35]); yticks(0:0.05:0.35); yticklabels(0:0.05*100:0.35*100)
xlabel('Pure tone CF (Hz)')
ylabel('Percent of samples')
legend({'Pitch Neurons','Non-Pitch Neurons'},'Location','northwest')
set(gca,'FontSize',20)

%%
% Low AAF (Ronnie 5)

PN_CFs = [];
nPN_CFs = [];

recording_idx = find(strcmp(PN_units(:,1),'Ronnie') & strcmp(PN_units(:,2),'P05'));
pen = recording_idx;
load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' PN_units{recording_idx,1}...
    '/tmp/Spikes_' PN_units{recording_idx,1} '_' PN_units{recording_idx,2} '_Good_Pitch.mat']);

stims = unique(type);
Flist = unique(F0);
repeats = unique(Y(:,5));
allUnits = unique(Y(:,3));
PNUnits = PN_units{pen,3};

[~,PNUnit_IDXs] = ismember(PNUnits,allUnits);
PN_tone_CFs = BFs(PNUnit_IDXs,13);
PN_CFs = [PN_CFs; PN_tone_CFs];
nPN_tone_CFs = BFs(:,13);
nPN_tone_CFs(PNUnit_IDXs) = [];
nPN_tone_CFs(nPN_tone_CFs==0) = [];
nPN_CFs = [nPN_CFs; nPN_tone_CFs];


figure;
[PN_N,PN_edges] = histcounts(PN_CFs,1:18);
PN_CFs_norm = PN_N/length(PN_CFs);

[nPN_N,~] = histcounts(nPN_CFs,1:18);
nPN_CFs_norm = nPN_N/length(nPN_CFs);
colors = colormap(jet(2));
histogram('Categories',string(Flist),'BinCounts',PN_CFs_norm,'FaceColor',colors(1,:),'FaceAlpha',0.5)
hold on
histogram('Categories',string(Flist),'BinCounts',nPN_CFs_norm,'FaceColor',colors(2,:),'FaceAlpha',0.5)
% ylim([0 0.30]); yticks(0:0.05:0.3); yticklabels(0:0.05*100:0.3*100)
xlabel('Pure tone CF (Hz)')
ylabel('Percent of samples')
legend({'Pitch Neurons','Non-Pitch Neurons'},'Location','northwest')
set(gca,'FontSize',20)

%%
% High AAF (Noah 8)

PN_CFs = [];
nPN_CFs = [];

recording_idx = find(strcmp(PN_units(:,1),'Noah') & strcmp(PN_units(:,2),'P08'));
pen = recording_idx;
load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' PN_units{recording_idx,1}...
    '/tmp/Spikes_' PN_units{recording_idx,1} '_' PN_units{recording_idx,2} '_Good_Pitch.mat']);

stims = unique(type);
Flist = unique(F0);
repeats = unique(Y(:,5));
allUnits = unique(Y(:,3));
PNUnits = PN_units{pen,3};

[~,PNUnit_IDXs] = ismember(PNUnits,allUnits);
PN_tone_CFs = BFs(PNUnit_IDXs,13);
PN_CFs = [PN_CFs; PN_tone_CFs];
nPN_tone_CFs = BFs(:,13);
nPN_tone_CFs(PNUnit_IDXs) = [];
nPN_tone_CFs(nPN_tone_CFs==0) = [];
nPN_CFs = [nPN_CFs; nPN_tone_CFs];


figure;
[PN_N,PN_edges] = histcounts(PN_CFs,1:18);
PN_CFs_norm = PN_N/length(PN_CFs);

[nPN_N,~] = histcounts(nPN_CFs,1:18);
nPN_CFs_norm = nPN_N/length(nPN_CFs);
colors = colormap(jet(2));
histogram('Categories',string(Flist),'BinCounts',PN_CFs_norm,'FaceColor',colors(1,:),'FaceAlpha',0.5)
hold on
histogram('Categories',string(Flist),'BinCounts',nPN_CFs_norm,'FaceColor',colors(2,:),'FaceAlpha',0.5)
ylim([0 0.60]); yticks(0:0.05:0.6); yticklabels(0:0.05*100:0.6*100)
xlabel('Pure tone CF (Hz)')
ylabel('Percent of samples')
legend({'Pitch Neurons','Non-Pitch Neurons'},'Location','northwest')
set(gca,'FontSize',20)




%%
% High A1 (Noah 3)

PN_CFs = [];
nPN_CFs = [];

recording_idx = find(strcmp(PN_units(:,1),'Noah') & strcmp(PN_units(:,2),'P03'));
pen = recording_idx;
load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' PN_units{recording_idx,1}...
    '/tmp/Spikes_' PN_units{recording_idx,1} '_' PN_units{recording_idx,2} '_Good_Pitch.mat']);

stims = unique(type);
Flist = unique(F0);
repeats = unique(Y(:,5));
allUnits = unique(Y(:,3));
PNUnits = PN_units{pen,3};

[~,PNUnit_IDXs] = ismember(PNUnits,allUnits);
PN_tone_CFs = BFs(PNUnit_IDXs,13);
PN_CFs = [PN_CFs; PN_tone_CFs];
nPN_tone_CFs = BFs(:,13);
nPN_tone_CFs(PNUnit_IDXs) = [];
nPN_tone_CFs(nPN_tone_CFs==0) = [];
nPN_CFs = [nPN_CFs; nPN_tone_CFs];


figure;
[PN_N,PN_edges] = histcounts(PN_CFs,1:18);
PN_CFs_norm = PN_N/length(PN_CFs);

[nPN_N,~] = histcounts(nPN_CFs,1:18);
nPN_CFs_norm = nPN_N/length(nPN_CFs);
colors = colormap(jet(2));
histogram('Categories',string(Flist),'BinCounts',PN_CFs_norm,'FaceColor',colors(1,:),'FaceAlpha',0.5)
hold on
histogram('Categories',string(Flist),'BinCounts',nPN_CFs_norm,'FaceColor',colors(2,:),'FaceAlpha',0.5)
ylim([0 0.60]); yticks(0:0.05:0.6); yticklabels(0:0.05*100:0.6*100)
xlabel('Pure tone CF (Hz)')
ylabel('Percent of samples')
legend({'Pitch Neurons','Non-Pitch Neurons'},'Location','northwest')
set(gca,'FontSize',20)

%% Mid A1 (Ronnie 13)

PN_CFs = [];
nPN_CFs = [];

recording_idx = find(strcmp(PN_units(:,1),'Ronnie') & strcmp(PN_units(:,2),'P13'));
pen = recording_idx;
load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' PN_units{recording_idx,1}...
    '/tmp/Spikes_' PN_units{recording_idx,1} '_' PN_units{recording_idx,2} '_Good_Pitch.mat']);

stims = unique(type);
Flist = unique(F0);
repeats = unique(Y(:,5));
allUnits = unique(Y(:,3));
PNUnits = PN_units{pen,3};

[~,PNUnit_IDXs] = ismember(PNUnits,allUnits);
PN_tone_CFs = BFs(PNUnit_IDXs,13);
PN_CFs = [PN_CFs; PN_tone_CFs];
nPN_tone_CFs = BFs(:,13);
nPN_tone_CFs(PNUnit_IDXs) = [];
nPN_tone_CFs(nPN_tone_CFs==0) = [];
nPN_CFs = [nPN_CFs; nPN_tone_CFs];


figure;
[PN_N,PN_edges] = histcounts(PN_CFs,1:18);
PN_CFs_norm = PN_N/length(PN_CFs);

[nPN_N,~] = histcounts(nPN_CFs,1:18);
nPN_CFs_norm = nPN_N/length(nPN_CFs);
colors = colormap(jet(2));
histogram('Categories',string(Flist),'BinCounts',PN_CFs_norm,'FaceColor',colors(1,:),'FaceAlpha',0.5)
hold on
histogram('Categories',string(Flist),'BinCounts',nPN_CFs_norm,'FaceColor',colors(2,:),'FaceAlpha',0.5)
ylim([0 0.35]); yticks(0:0.05:0.35); yticklabels(0:0.05*100:0.35*100)
xlabel('Pure tone CF (Hz)')
ylabel('Percent of samples')
legend({'Pitch Neurons','Non-Pitch Neurons'},'Location','northwest')
set(gca,'FontSize',20)


%
% histogram(PN_CFs,1:18,'FaceColor',colors(1,:),'FaceAlpha',0.5)
% hold on
% histogram(nPN_CFs,1:18,'FaceColor',colors(2,:),'FaceAlpha',0.5)
% xticks(1:18)
% xticklabels(Flist)
% xlabel('Pure tone BF (Hz')
% ylabel('Number of neurons')
% legend({'Pitch Neurons','Non-Pitch Neurons'},'Location','northwest')
% set(gca,'FontSize',20)
% title('High AAF (Noah 8)')
% title('High A1 (Noah 3)')
% title('Mid A1 (Ronnie 13)')
