%% Function to plot tuning curves of specified stimuli for a single unit
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, November 2022

function a = plot_tuning_by_cond(Y,type,F0,unit,stims,BFs,animal,pen)

    % INPUTS:
    % Y - n x 6 matrix where Y(:,)=spiketimes, Y(:,3)=unit spiking, Y(:,4)=stimulus, Y(:,5)=repeat, Y(:,6)=trial number
    % type - nStimuli x 1 vector labeling the sound type (CT0, tone, etc) of each stimulus
    % F0 - nStimuli x 1 vector labeling the F0 of each stimulus
    % unit - the unit to be plotted
    % stims - cell array of strings containing stimuli names to be plotted
    % BFs - length(stims) x 1 with the input unit's best frequencies for each of the sound types listed in 'stims'
    % animal - string, animal this unit came from
    % pen - string, penetration this unit came from

    Flist = unique(F0);
    repeats = unique(Y(:,5));
    window = [0 0.8];
    colors = colormap(hsv(length(stims))); % make the colormap to be used later
    
    figure('Position',[1900 500 1800 1200])

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
        meanSpikes = mean(nSpikes); % average across repeats

        %%%%%%%%%%%%%%%%%%%%% plot the PTSH %%%%%%%%%%%%%%%%%%%%%%%%%

        h = shadedErrorBar(1:17,meanSpikes,ste(nSpikes),{'Color',colors(ss,:)},1);
        h.mainLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
        h.edge(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
        h.edge(2).Annotation.LegendInformation.IconDisplayStyle = 'off';
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
        set(gca,'fontsize',22)
        axis tight

    end

    legend(stims)
    
    return
end