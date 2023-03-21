%% Plot all the responses to all stimuli for one unit

% Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
%     'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
%     'Dory','Dory','Dory','Dory'};
% 
% Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
%     'P04','P05','P08','P13','P02','P03','P05','P08',...
%     'P00','P01','P02','P04'};
% 
% Qualia = 'Good';

load('/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/Derry/Spikes_Derry_P05_Good_Pitch.mat');

units = unique(Y(:,3));

for uu = 4:10
    unit = units(uu);

    thisUnitY = Y(Y(:,3)==unit,:);

    stims = unique(thisUnitY(:,4));

    rasterHeightCounter = 0;

    for ss = 1:length(stims)

        stim = stims(ss);
        repeats = unique(thisUnitY(thisUnitY(:,4)==stim,5));
        
        for rr = 1:length(repeats)

            repeat = repeats(rr);
            trialRows = thisUnitY(:,4)==stim & thisUnitY(:,5)==repeat;
            scatter(thisUnitY(trialRows,2),ones(length(find(trialRows)),1)*rasterHeightCounter,'k.');
            hold on
    
            rasterHeightCounter = rasterHeightCounter+1;

        end
    
    end

    axis tight
    hold off
    pause

end


%% Plot the range of values in each column of Y

units = unique(Y(:,3));
maxStim = zeros(length(units),1);
maxRepeat = zeros(length(units),1);
maxSweep = zeros(length(units),1);

for uu = 1:length(units)

    thisUnit = Y(Y(:,3)==units(uu),:);

    maxStim(uu) =  max(thisUnit(:,4));
    maxRepeat(uu) = max(thisUnit(:,5));
    maxSweep(uu) = max(thisUnit(:,6));

end

figure; hist(maxStim);
figure; hist(maxRepeat);
figure; hist(maxSweep);


%% Plot response profiles for one unit at a time, to many stims

load('/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/Noah/Spikes_Noah_P02_Good_Pitch.mat');

stims = {'CT0','CT5','CT20','CT40'};
Flist = unique(F0);
repeats = unique(Y(:,5));
units = unique(Y(:,3));
window = [0 0.5];

figure('Position',[1500 500 1800 1200])
for uu = 1:length(units)

    clf;
    sgtitle(sprintf('Unit # %d',units(uu)))

    unitSpikes = Y(Y(:,3)==181,:);
    
    % go through each stim we want to plot
    for ss = 1:length(stims)

        nSpikes = zeros(length(repeats),length(Flist)); %creates an array of length repList by length Flist
        spikeTimes = cell(1,length(Flist)); %blank array 
            
        for ff = 1:length(Flist) % do the below for all the frequencies
    
            stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); 

            if isempty(stimNum) % if this stim type and fo combo wasn't presented
                nSpikes(ff) = 0;
                spikes = 0;
                continue
            end
            
            % finds the stim label that corresponds to this stim type at
            % this particular F0

            % this stimNum will have been presented multiple times, 
            % so go through each presentation
    
            spikes = [];
            for rr = 1:length(repeats)
            
                spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                nSpikes(rr,ff) = sum(spikeIDXs);
                spikes = [spikes; unitSpikes(spikeIDXs,2)];

            end

            spikeTimes{ff} = spikes;
            
        end

        nSpikes = nSpikes ./ diff(window); % spikes per second
        spikes = sum(nSpikes,'all');
        meanspikes = spikes/(length(Flist)*length(repeats));

        %%%%%%%%%%%%%%%%%%% plot the PTSH %%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(2,4,(ss*2));
     
        errorbar(1:length(Flist),mean(nSpikes),ste(nSpikes)) 
        % error bar plot for X axis frequency, y axis mean spike number, error bars are standard error for the mean spike number
    
        xticks([])
        ylabel('Evoked Firing Rate (spike/sec)'); % y axis label is the firing rate
        set(gca,'XDir','reverse')
        camroll(-90)
        axis tight
      
        %%%%%%%%%%%%%%%%% plot the rasters %%%%%%%%%%%%%%%%%%%%%%%
        subplot(2,4,(ss*2)-1); hold on
        title(stims(ss))

        for k = 1:length(spikeTimes) 
            s = scatter(spikeTimes{k},ones(length(spikeTimes{k}),1)*k,'marker','o','markerfacecolor','k','markeredgecolor','k'); %create a scatter plot of spike times versus
            s.SizeData = 10; %make the dots in the raster plot that size
        end

        yticks(1:length(Flist))
        yticklabels(string(Flist))
        
        hold off
        xlabel('Time relative to stimulus onset (second)')
        ylabel('F0 (Hz)')
        ylim([1 length(Flist)])

    end
    pause
end


%% Display the number of good units per penetration 

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';



for rr = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{rr} '/Spikes_' Animals{rr} '_' Pens{rr} '_Good_Pitch.mat']);


    sprintf('%s, %s: %d neurons', Animals{rr}, Pens{rr}, length(unique(Y(:,3))))
end


%% Display the ISI histogram and latency by F0 across multiple conditions

load('/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/Derry/Spikes_Derry_P05_Good_Pitch.mat');

stims = {'CT0','low','high','tone'};
Flist = unique(F0);
repeats = unique(Y(:,5));
units = unique(Y(:,3));
window = [0 0.5];

figure('Position',[1500 500 1800 1200])
for uu = 1:length(units)

    clf;
    sgtitle(sprintf('Unit # %d',units(uu)))

    unitSpikes = Y(Y(:,3)==units(uu),:);

    % go through each stim we want to plot
    for ss = 1:length(stims)
            
        ISIs = cell(length(Flist),1);
        onsetDelays = zeros(length(repeats), length(Flist));

        for ff = 1:length(Flist) % do the below for all the frequencies
    
            stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); 

            if isempty(stimNum) % if this stim type and fo combo wasn't presented
                 continue
            end
            
            % finds the stim label that corresponds to this stim type at
            % this particular F0

            % this stimNum will have been presented multiple times, 
            % so go through each presentation
    
            spikeTimeDiffs = [];
            for rr = 1:length(repeats)
                % take all the spike times of this particular stimulation
                spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);

                % get the onset of the first spike and add it to the vector
                if ~isempty(find(spikeIDXs))

                    relativeSpikes = unitSpikes(spikeIDXs,2);
                    ISIvector = diff(relativeSpikes);
                    spikeTimeDiffs = [spikeTimeDiffs; ISIvector];
        
                    onsetDelays(rr,ff) = relativeSpikes(1);
                end       
            end
     
            ISIs{ff} = spikeTimeDiffs;

        end

        %%%%%%%%%%%%%%%%%%% plot the ISI histogram %%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(2,4,(ss*2));

        colors = colormap(parula(17));
        for ff=1:length(Flist)

            histogram(ISIs{ff},'FaceColor',colors(ff,:),'FaceAlpha',0.2,'BinEdges',0.005:0.005:0.1)
            hold on

        end

        xlabel('ISI')
    
      
        %%%%%%%%%%%%%%%%% plot the latency histogram %%%%%%%%%%%%%%%%%%%%%%%
        subplot(2,4,(ss*2)-1); hold on
        title(stims(ss))

        for ff=1:length(Flist)

            histogram(onsetDelays(:,ff),'FaceColor',colors(ff,:),'FaceAlpha',0.2,'BinEdges',0.005:0.005:0.1)
            hold on

        end

        xlabel('Spiking onset relative to stimulus')
        
    end
    
    pause

end


%% Display the onset and ISI histograms for each penetration

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

tiledlayout(5,4);

for rr = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{rr} '/Spikes_' Animals{rr} '_' Pens{rr} '_Good_Pitch.mat']);

    sweeps = unique(Y(:,6));
    onsets = [];
    ISIs = [];

    for ss = 1:length(sweeps)

        thisSweep = Y(Y(:,6)==sweeps(ss),2);

        onsets = [onsets; thisSweep(1)];
        ISIs = [ISIs; diff(thisSweep)];

    end

    nexttile;
    histogram(onsets,'BinEdges',0:0.01:0.3);
    
    title([Animals(rr) ' ' Pens(rr)])
end
