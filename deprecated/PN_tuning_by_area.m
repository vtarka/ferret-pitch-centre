%% Veronica Tarka
% veronica.tarka@dpag.ox.ac.uk
% January 2023

%%
% High A1 (Noah 3)

recording_idx = find(strcmp(PN_units(:,1),'Noah') & strcmp(PN_units(:,2),'P08'));
pen = recording_idx;
load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' PN_units{recording_idx,1}...
    '/tmp/Spikes_' PN_units{recording_idx,1} '_' PN_units{recording_idx,2} '_Good_Pitch.mat']);

Flist = unique(F0);
repeats = unique(Y(:,5));
window = [0 0.15];
allUnits = unique(Y(:,3));
PNUnits = PN_units{pen,3};
[~,PNUnit_IDXs] = ismember(PNUnits,allUnits);

nPNUnits = allUnits;
nPNUnits(PNUnit_IDXs) = [];

nPN_tone_CFs = BFs(:,13);
nPN_tone_CFs(PNUnit_IDXs) = [];
nPNUnits(nPN_tone_CFs==0) = [];
[~,nPNUnit_IDXs] = ismember(nPNUnits,allUnits);


stims_to_plot = {'tone'};
colors = colormap(jet(2));

figure('Position',[1500 500 1800 1200])

for uu = 1:length(PNUnits)
    unit = PNUnits(uu);

    unitSpikes = Y(Y(:,3)==unit,:); % get the spikes of just this unit

    tuning = zeros(length(stims_to_plot),length(Flist));

    for ss = 1:length(stims_to_plot)
        
        nSpikes = zeros(length(repeats),length(Flist)); %creates an array of length repList by length Flist
        
        for ff = 1:length(Flist) % do the below for all the frequencies
    
            stimNum = find(strcmp(type,stims_to_plot{ss}) & (F0==Flist(ff))); 

            if isempty(stimNum) % if this stim type and fo combo wasn't presented
                nSpikes(:,ff) = 0;
                continue
            end
            
            % finds the stim label that corresponds to this stim type at
            % this particular F0

            % this stimNum will have been presented multiple times, 
            % so go through each presentation
    
            for rr = 1:length(repeats)
            
                spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                nSpikes(rr,ff) = sum(spikeIDXs);

            end   
        end % ends F0 loop

        nSpikes = nSpikes ./ diff(window); % spikes per second

        meanSpikes = mean(nSpikes);

        subplot(2,4,uu)
        h = shadedErrorBar(1:17,meanSpikes,ste(nSpikes),{'Color',colors(ss,:)},1);
        h.mainLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
        h.edge(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
        h.edge(2).Annotation.LegendInformation.IconDisplayStyle = 'off';
        hold on
        if BFs(ss)~= 0
            p = scatter(BFs(ss),meanSpikes(BFs(ss)),200,'MarkerFaceColor',colors(ss,:));
            p.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end


        xticks(1:17)
        xticklabels(num2str(Flist))
        axis tight
        xlabel('F0')
        ylabel('Evoked Firing Rate (spike/sec)'); % y axis label is the firing rate

        title(sprintf('BF=%d Hz',Flist(BFs(PNUnit_IDXs(uu),13))))


    end % end stim loop
end

%%
% High A1 (Noah 3)

recording_idx = find(strcmp(PN_units(:,1),'Noah') & strcmp(PN_units(:,2),'P02'));
pen = recording_idx;
load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' PN_units{recording_idx,1}...
    '/tmp/Spikes_' PN_units{recording_idx,1} '_' PN_units{recording_idx,2} '_Good_Pitch.mat']);

stims = unique(type);
Flist = unique(F0);
repeats = unique(Y(:,5));
window = [0 0.15];
allUnits = unique(Y(:,3));
PNUnits = PN_units{pen,3};
[~,PNUnit_IDXs] = ismember(PNUnits,allUnits);

nPNUnits = allUnits;
nPNUnits(PNUnit_IDXs) = [];

nPN_tone_CFs = BFs(:,13);
nPN_tone_CFs(PNUnit_IDXs) = [];
nPNUnits(nPN_tone_CFs==0) = [];

stims_to_plot = {'tone'};
colors = colormap(jet(2));

figure; 
for uu = 1:length(allUnits)
    unit = allUnits(uu);

    unitSpikes = Y(Y(:,3)==unit,:); % get the spikes of just this unit

    tuning = zeros(length(stims_to_plot),length(Flist));

    for ss = 1:length(stims_to_plot)
        
        nSpikes = zeros(length(repeats),length(Flist)); %creates an array of length repList by length Flist
        
        for ff = 1:length(Flist) % do the below for all the frequencies
    
            stimNum = find(strcmp(type,stims_to_plot{ss}) & (F0==Flist(ff))); 

            if isempty(stimNum) % if this stim type and fo combo wasn't presented
                nSpikes(:,ff) = 0;
                continue
            end
            
            % finds the stim label that corresponds to this stim type at
            % this particular F0

            % this stimNum will have been presented multiple times, 
            % so go through each presentation
    
            for rr = 1:length(repeats)
            
                spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                nSpikes(rr,ff) = sum(spikeIDXs);

            end   
        end % ends F0 loop

        nSpikes = nSpikes ./ diff(window); % spikes per second

        meanSpikes = mean(nSpikes);

        figure; 
        h = shadedErrorBar(1:17,meanSpikes,ste(nSpikes),{'Color',colors(ss,:)},1);
        h.mainLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
        h.edge(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
        h.edge(2).Annotation.LegendInformation.IconDisplayStyle = 'off';
        hold on
        if BFs(ss)~= 0
            p = scatter(BFs(ss),meanSpikes(BFs(ss)),200,'MarkerFaceColor',colors(ss,:));
            p.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end


        xticks(1:17)
        xticklabels(num2str(Flist))
        xlabel('F0')
        ylabel('Evoked Firing Rate (spike/sec)'); % y axis label is the firing rate

        title(sprintf('BF=%d',BFs(uu,13)))
        pause
        close all

    end % end stim loop
end