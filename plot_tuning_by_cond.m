%% Function to plot tuning curves of specified stimuli for a single unit
% DEPENDENCIES: shadedErrorBar.m
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, November 2022

function ymax = plot_tuning_by_cond(Y,type,F0,unit,stims,BFs,animal,pen,window,colors)

    % INPUTS:
    % Y - n x 6 matrix where Y(:,)=spiketimes, Y(:,3)=unit spiking, Y(:,4)=stimulus, Y(:,5)=repeat, Y(:,6)=trial number
    % type - nStimuli x 1 vector labeling the sound type (CT0, tone, etc) of each stimulus
    % F0 - nStimuli x 1 vector labeling the F0 of each stimulus
    % unit - the unit to be plotted
    % stims - cell array of strings containing stimuli names to be plotted
    % BFs - length(stims) x 1 with the input unit's best frequencies for each of the sound types listed in 'stims'
    % animal - string, animal this unit came from
    % pen - string, penetration this unit came from
    % window - vector with two elements: window(1) is the start of the response window, window(2) is the end (in ms) (OPTIONAL)
    % colors - nStims x 3 matrix, each row holds an RGB triplet of the color for that stimulus (OPTIONAL)

    Flist = unique(F0);
    repeats = unique(Y(:,5));
    if ~exist('window','var')
        window = [0 0.08];
    end

    if ~exist('colors','var')
        colors = colormap(hsv(length(stims))); % make the colormap to be used later
    end
    
%     figure('Position',[1900 500 1800 1200])
%     figure
    sgtitle(sprintf('%s, %s unit # %d',animal,pen,unit)) % label the plot with the animal, penetration, and unit shown

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
        nSpikes = zscore(nSpikes,0,'all');
        meanSpikes = zscore(mean(nSpikes)); % average across repeats

%         spike_error = ste(nSpikes);

        %%%%%%%%%%%%%%%%%%%%% plot the PTSH %%%%%%%%%%%%%%%%%%%%%%%%%

        if strcmp('low',stims{ss}) || strcmp('F0MaskLow',stims{ss})
            if contains(stims{ss},'Mask')
                h = shadedErrorBar(3:17,meanSpikes(3:end),ste(nSpikes(:,3:end)),{'Color',colors(ss,:),'linestyle','--'},1);
            else
                h = shadedErrorBar(3:17,meanSpikes(3:end),ste(nSpikes(:,3:end)),{'Color',colors(ss,:),'linewidth',3},1);
            end
        elseif contains(stims{ss},'tone')
            h = shadedErrorBar(1:17,meanSpikes,ste(nSpikes),{'Color',colors(ss,:),'linestyle','--','linewidth',3},1);
        else
            h = shadedErrorBar(1:17,meanSpikes,ste(nSpikes),{'Color',colors(ss,:),'linewidth',3},1);
        end

%         if length(stims)==2 && ss == 2
%             h = shadedErrorBar(1:17,meanSpikes,ste(nSpikes),{'Color',colors(3,:)},1);
%         else
%             h = shadedErrorBar(1:17,meanSpikes,ste(nSpikes),{'Color',colors(ss,:)},1);
%         end

%         h = shadedErrorBar(1:17,meanSpikes,ste(nSpikes),{'Color',colors(ss,:)},1);
        h.mainLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
        h.edge(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
        h.edge(2).Annotation.LegendInformation.IconDisplayStyle = 'off';
        hold on
%         errorbar(1:17,meanSpikes,ste(nSpikes),'color',colors(ss,:));
        hold on
        
        % if this unit is frequency sensitive to this sound type
        if BFs(ss)~= 0
            p = scatter(BFs(ss),meanSpikes(BFs(ss)),200,'MarkerFaceColor',colors(ss,:)); % add the best frequency for this sound type as a point
            p.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end

        xticks(1:17)
        xticklabels(num2str(Flist))
        xlabel('F0')
        ylabel('Evoked Firing Rate (spike/sec)'); % y axis label is the firing rate
%         set(gca,'fontsize',22)
        axis tight

    end

%     legend(stims)
    
    return
end