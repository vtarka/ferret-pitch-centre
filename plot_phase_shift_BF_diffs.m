%% Plot histogram of tuning differences between high and alt, and high and rand
% DEPENDENCIES: plot_tuning_by_cond.m
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, March 2023

load('HNs_fNoah.mat')
load('TNs_fNoah.mat')
load('PNs_fNoah.mat')

HNalt_diffs = [];
HNrand_diffs = [];

for pen = 1:size(HN_units,1)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' HN_units{pen,1} '/tmp02/Spikes_' HN_units{pen,1} '_' HN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'high','alt','rand'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));

    HNUs = HN_units{pen,3};

    for uu = 1:size(HNUs,1)

        if HNUs(uu,2)==1
            window = [0 0.06];
        elseif HNUs(uu,2)==2
            window = [0.06 0.15];
        else
            if pen<9
                window = [.3 .4];
            else
                window = [.2 .3];
            end
        end

        BFs = zeros(3,1);

        unitSpikes = Y(Y(:,3)==HNUs(uu,1),:); % extract the spikes of just this unit

        % for each stim type
        for ss = 1:length(stims)
                
            meanSpikes = zeros(length(Flist),1);

            % go through each F0 to find BF
            for ff = 1:length(Flist)

                stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); % unique name for combination of stim type and F0

                if isempty(stimNum) % if this stim type and fo combo wasn't presented
                    continue
                end

                nSpikes = zeros(length(repeats),1); % initialize space to count the number of spikes per trial repeat

                % for each repeat
                for rr = 1:length(repeats)

                    spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                    nSpikes(rr) = sum(spikeIDXs);

                end % ends repeat loop

                meanSpikes(ff) = mean(nSpikes); % average the spike rates across repeats

            end % ends looping through F0s
            
            [~,i] = max(meanSpikes); % find the index where the maximum spike rate occurs
            BFs(ss) = i; % save the index of the best frequency
        end

        HNalt_diffs(end+1) = BFs(1) - BFs(2);
        HNrand_diffs(end+1) = BFs(1) - BFs(3);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%% TNs %%%%%%%%%%%%%%%%%%%%%%%%%%%

TNalt_diffs = [];
TNrand_diffs = [];

for pen = 1:size(TN_units,1)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' TN_units{pen,1} '/tmp02/Spikes_' TN_units{pen,1} '_' TN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'high','alt','rand'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));

    TNUs = TN_units{pen,3};

    for uu = 1:size(TNUs,1)

        if TNUs(uu,2)==1
            window = [0 0.06];
        elseif TNUs(uu,2)==2
            window = [0.06 0.15];
        else
            if pen<9
                window = [.3 .4];
            else
                window = [.2 .3];
            end
        end

        BFs = zeros(3,1);

        unitSpikes = Y(Y(:,3)==TNUs(uu,1),:); % extract the spikes of just this unit

        % for each stim type
        for ss = 1:length(stims)
                
            meanSpikes = zeros(length(Flist),1);

            % go through each F0 to find BF
            for ff = 1:length(Flist)

                stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); % unique name for combination of stim type and F0

                if isempty(stimNum) % if this stim type and fo combo wasn't presented
                    continue
                end

                nSpikes = zeros(length(repeats),1); % initialize space to count the number of spikes per trial repeat

                % for each repeat
                for rr = 1:length(repeats)

                    spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                    nSpikes(rr) = sum(spikeIDXs);

                end % ends repeat loop

                meanSpikes(ff) = mean(nSpikes); % average the spike rates across repeats

            end % ends looping through F0s
            
            [~,i] = max(meanSpikes); % find the index where the maximum spike rate occurs
            BFs(ss) = i; % save the index of the best frequency
        end

        TNalt_diffs(end+1) = BFs(1) - BFs(2);
        TNrand_diffs(end+1) = BFs(1) - BFs(3);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PNs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PNalt_diffs = [];
PNrand_diffs = [];

for pen = 1:size(PN_units,1)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' PN_units{pen,1} '/tmp02/Spikes_' PN_units{pen,1} '_' PN_units{pen,2} '_Good_Pitch.mat']);

    stims = {'high','alt','rand'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));

    PNUs = PN_units{pen,3};

    for uu = 1:size(PNUs,1)

        if PNUs(uu,2)==1
            window = [0 0.06];
        elseif PNUs(uu,2)==2
            window = [0.06 0.15];
        else
            if pen<9
                window = [.3 .4];
            else
                window = [.2 .3];
            end
        end

        BFs = zeros(3,1);

        unitSpikes = Y(Y(:,3)==PNUs(uu,1),:); % extract the spikes of just this unit

        % for each stim type
        for ss = 1:length(stims)
                
            meanSpikes = zeros(length(Flist),1);

            % go through each F0 to find BF
            for ff = 1:length(Flist)

                stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); % unique name for combination of stim type and F0

                if isempty(stimNum) % if this stim type and fo combo wasn't presented
                    continue
                end

                nSpikes = zeros(length(repeats),1); % initialize space to count the number of spikes per trial repeat

                % for each repeat
                for rr = 1:length(repeats)

                    spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                    nSpikes(rr) = sum(spikeIDXs);

                end % ends repeat loop

                meanSpikes(ff) = mean(nSpikes); % average the spike rates across repeats

            end % ends looping through F0s
            
            [~,i] = max(meanSpikes); % find the index where the maximum spike rate occurs
            BFs(ss) = i; % save the index of the best frequency
        end

        PNalt_diffs(end+1) = BFs(1) - BFs(2);
        PNrand_diffs(end+1) = BFs(1) - BFs(3);
    end
end



%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%

figure;

subplot(1,3,1)
histogram(HNalt_diffs)
hold on
histogram(HNrand_diffs)
set(gca,'fontsize',24)
title('Harmonicity Neurons')


subplot(1,3,2)
histogram(TNalt_diffs)
hold on
histogram(TNrand_diffs)
set(gca,'fontsize',24)
title('Periodicity Neurons')


subplot(1,3,3)
histogram(PNalt_diffs)
hold on
histogram(PNrand_diffs)
set(gca,'fontsize',24)
legend('high BF - alt BF','high BF - rand BF')
title('Pitch Neurons')