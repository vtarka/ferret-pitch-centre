%% Plot the tuning curves evoked by stimulus onset AND offset if the neuron showed an MI peak following either
% DEPENDENCIES:
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

timebins = 0:15:500;
uCounter = 1;
colors = colormap(hsv(13));
figure('Position',[1900 500 1800 1200])

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    for uu = 1:length(units) % for each unit
        unit = units(uu);
        unitSpikes = Y(Y(:,3)==unit,:); % get the spikes of just this unit

        if isempty(find(all_peaks(uCounter,:,:), 1))
            uCounter = uCounter + 1;
            continue
        end

        clf;
        stims4onsetlegend = [];
        stims4offsetlegend = [];
        for ss = 1:length(stims)

            % if this stimulus evoked a peak in the MI timecourse within 150 ms, evaluate the onset tuning
            if all_peaks(uCounter, ss, 1) ~=0 && all_peaks(uCounter,ss,1) < 12 % the 12th timebin corresponds to 165 ms
                stims4onsetlegend = [stims4onsetlegend; ss];
                window = [0 timebins(all_peaks(uCounter,ss,1)+1)/1000];

                nSpikes = zeros(length(repeats),length(Flist)); % allocate space to save spiking info for each trial

                for ff = 1:length(Flist)

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

                end % ends F0 loop

                nSpikes = nSpikes ./ diff(window); % spikes per second
                meanSpikes = mean(nSpikes); % average across repeats

                subplot(2,1,1); hold on
                plot(1:17, meanSpikes,'Color',colors(ss,:),'LineWidth',3);
            end % ends onset if condition

            % if this stimulus evoked a peak in the MI timecourse within 150 ms after offset, evaluate the offset tuning
            if strcmp(Animals{ap},'Noah')
                if all_peaks(uCounter, ss, 2) ~=0 && all_peaks(uCounter,ss,2) < 32 % the 32nd timebin corresponds to 465 ms
                    stims4offsetlegend = [stims4offsetlegend; ss];
                    window = [0.3 timebins(all_peaks(uCounter,ss,2)+1)/1000];
    
                    nSpikes = zeros(length(repeats),length(Flist)); % allocate space to save spiking info for each trial
    
                    for ff = 1:length(Flist)
    
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
    
                    end % ends F0 loop
    
                    nSpikes = nSpikes ./ diff(window); % spikes per second
                    meanSpikes = mean(nSpikes); % average across repeats
    
                    subplot(2,1,2); hold on
                    plot(1:17, meanSpikes,'Color',colors(ss,:),'LineWidth',3);
                end
            else
                if all_peaks(uCounter, ss, 2) ~=0 && all_peaks(uCounter,ss,2) < 25 % the 25th timebin corresponds to 345 ms
                    stims4offsetlegend = [stims4offsetlegend; ss];
                    window = [0.2 timebins(all_peaks(uCounter,ss,2)+1)/1000];
    
                    nSpikes = zeros(length(repeats),length(Flist)); % allocate space to save spiking info for each trial
    
                    for ff = 1:length(Flist)
    
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
    
                    end % ends F0 loop
    
                    nSpikes = nSpikes ./ diff(window); % spikes per second
                    meanSpikes = mean(nSpikes); % average across repeats
    
                    subplot(2,1,2); hold on
                    plot(1:17, meanSpikes,'Color',colors(ss,:),'LineWidth',3);

                end % ends offset if condition

            end % ends Noah if case

        end % ends stimulus loop

        subplot(2,1,1)
        title('Onset')
        legend(stims(stims4onsetlegend))
        axis tight
        set(gca,'FontSize',22)

        subplot(2,1,2)
        title('Offset')
        legend(stims(stims4offsetlegend))
        axis tight
        set(gca,'FontSize',22)

        sgtitle(sprintf('%s %s Unit: %d',Animals{ap},Pens{ap},unit))

        pause

        uCounter = uCounter + 1;
    end % ends unit loop

end % ends recording loop