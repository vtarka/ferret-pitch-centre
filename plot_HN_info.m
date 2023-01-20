
%% Make plots on HNs including tuning curves, where they are, etc
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, January 2023

load('HN_units.mat')

%% First plot all of them

figure('Position',[1900 500 1800 1200])
stims = {'high','low','CT0'};
colors = colormap(hsv(length(stims))); % make the colormap to be used later
sp = 0;

% for each penetration
for pen = 1:length(HN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' HN_units{pen,1} '/tmp/Spikes_' HN_units{pen,1} '_' HN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'high','low','CT0'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    window = [0 0.15]; % response window

    HNUnits = HN_units{pen,3}; % array of units we found to be harmonicity neurons
    [~,HNUnit_IDXs] = ismember(HNUnits,allUnits); % find the indices of these pitch neurons within the list of all units

    HN_BFs = BFs(HNUnit_IDXs,13); % get the PN's best frequencies

    for hn = 1:length(HNUnits)

        sp = sp + 1;

        if sp > 10
            break
        end

        unit = HNUnits(hn);

        unitSpikes = Y(Y(:,3)==unit,:); % get spikes for just this unit
    
        % go through each stim we want to plot
        for ss = 1:length(stims)
    
            nSpikes = zeros(length(repeats),length(Flist)); % allocate space to save spiking info for each trial
                
            for ff = 1:length(Flist) % for each frequency
        
                stimNum = find(strcmp(type,stims{ss}) & (F0==Flist(ff))); % unique name for combination of stim type and F0
    
                if isempty(stimNum) % if this stim type and fo combo wasn't presented
                    nSpikes(:,ff) = 0;
                    continue
                end
        
                % for each trial of this stim type
                for rr = 1:length(repeats)
                
                    spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                    nSpikes(rr,ff) = sum(spikeIDXs);
    
                end   
            end
    
            nSpikes = nSpikes ./ diff(window); % spikes per second
            meanSpikes = mean(nSpikes); % average across repeats
    
            %%%%%%%%%%%%%%%%%%%%% plot the PTSH %%%%%%%%%%%%%%%%%%%%%%%%%
    
         
            subplot(5,2,sp)

            if ss ==  1
                plot(1:17,meanSpikes,'Color',colors(ss,:),'LineWidth',3,'LineStyle',':')
            else
                plot(1:17,meanSpikes,'Color',colors(ss,:),'LineWidth',3)
            end

            hold on
    
            axis tight

            if sp==9
                xticks(5:5:20)
                xticklabels(Flist(5:5:17))
                xlabel('F0 (Hz)')

                yticks(0:5)
                ylabel('Spks/s')
            else
                xticks([])
                yticks(0:10:60)

                xticks(5:5:20)
                xticklabels(Flist(5:5:17))
            end

            set(gca,'FontSize',20)


        end % ends stim loop

    end % ends unit loop
end % ends recording loop


%% Next evaluate locations

penUnits = [27,146,156,211,36,46,9,43,...
    34,77,146,58,60,61,37,29,...
    43,40,26,26]; % the numer of total units in each penetration (in each row of PN_units)

% Make a list of how many PNs per penetration there were
HNUnits = [];

for ap = 1:length(HN_units)
    HNUnits = [HNUnits; length(HN_units{ap,3})];
end

figure;
HNUnitFrequency = HNUnits; %./ sum(HNUnits); % percent of units in each penetration determined to be PN_units
plot(HNUnitFrequency,'Color','k','LineWidth',3)
xticks(1:20);xticklabels({'N1','N2','N3','N4','N5','N6','N7','N8','R4','R5','R8','R13','De2','De3','De5','De8','Do0','Do1','Do2','Do4'})
% ylim([0 0.35]); yticks(0:0.05:0.35); yticklabels(0:0.05*100:0.35*100)
ylabel('Percent of neurons found to be HNs')
xlabel('Penetration')
set(gca,'FontSize',18)
xlim([1 20])


%% Next evaluate the differences in best frequency differences between CT0 and low harmonics

BF_diffs = [];

% for each penetration
for pen = 1:length(HN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' HN_units{pen,1} '/tmp/Spikes_' HN_units{pen,1} '_' HN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'high','low','CT0'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    window = [0 0.15]; % response window

    HNUnits = HN_units{pen,3}; % array of units we found to be harmonicity neurons
    [~,HNUnit_IDXs] = ismember(HNUnits,allUnits); % find the indices of these pitch neurons within the list of all units

    HN_BFs = BFs(HNUnit_IDXs,:); % get the PN's best frequencies

    for hn = 1:length(HNUnits)
        BF_diffs = [BF_diffs; abs(HN_BFs(hn,1)-HN_BFs(hn,11))];
    end

end

figure; histogram(BF_diffs)


%% Evaluate the peakedness

figure('Position',[1900 500 1800 1200])
stims = {'CT0'};
colors = colormap(hsv(length(stims))); % make the colormap to be used later

peak_counts = [];
peak_diffs = [];

% for each penetration
for pen = 1:length(HN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' HN_units{pen,1} '/tmp/Spikes_' HN_units{pen,1} '_' HN_units{pen,2} '_Good_Pitch.mat']);

    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    window = [0 0.15]; % response window

    HNUnits = HN_units{pen,3}; % array of units we found to be harmonicity neurons
    [~,HNUnit_IDXs] = ismember(HNUnits,allUnits); % find the indices of these pitch neurons within the list of all units

    for hn = 1:length(HNUnits)

        unit = HNUnits(hn);

        unitSpikes = Y(Y(:,3)==unit,:); % get spikes for just this unit
    
        % go through each stim we want to plot
        for ss = 1:length(stims)
    
            nSpikes = zeros(length(repeats),length(Flist)); % allocate space to save spiking info for each trial
                
            for ff = 1:length(Flist) % for each frequency
        
                stimNum = find(strcmp(type,stims{ss}) & (F0==Flist(ff))); % unique name for combination of stim type and F0
    
                if isempty(stimNum) % if this stim type and fo combo wasn't presented
                    nSpikes(:,ff) = 0;
                    continue
                end
        
                % for each trial of this stim type
                for rr = 1:length(repeats)
                
                    spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                    nSpikes(rr,ff) = sum(spikeIDXs);
    
                end   
            end
    
            clf
            nSpikes = nSpikes ./ diff(window); % spikes per second
            meanSpikes = mean(nSpikes); % average across repeats

            [nPeaks, peak_idx] = count_peaks(meanSpikes,0.75);

            peak_counts = [peak_counts ; nPeaks];
            peak_diffs = [peak_diffs; diff(peak_idx)];

%             plot(1:17,meanSpikes,'LineWidth',3)
%             hold on
% 
%             for peak = 1:length(peak_idx)
%                 scatter(peak_idx(peak),meanSpikes(peak_idx(peak)),200,'filled')
%             end
% 
%             pause

        end
    end
end

