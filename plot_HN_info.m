
%% Make plots on HNs including tuning curves, where they are, etc
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, January 2023

load('HN_units.mat')

%% First plot all of them

figure('Position',[1900 500 1800 1200])
stims = {'high','low','CT0','tone'};
colors = colormap(hsv(length(stims))); % make the colormap to be used later
sp = 0;

% for each penetration
for pen = 1:length(HN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' HN_units{pen,1} '/tmp/Spikes_' HN_units{pen,1} '_' HN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'high','low','CT0','tone'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    HNUnits = HN_units{pen,3}; % array of units we found to be harmonicity neurons
    [~,HNUnit_IDXs] = ismember(HNUnits,allUnits); % find the indices of these pitch neurons within the list of all units

    HN_BFs = BFs(HNUnit_IDXs,13); % get the PN's best frequencies

    for hn = 1:length(HNUnits)

        sp = sp + 1;

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

