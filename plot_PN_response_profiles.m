%% Plot a 2D matrix with each row a z-scored response to a given stimulus
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, March 2023

load('PNs_fNoah.mat')

stims = {'low','high','CT0','allHarm','tone'};
windows = [0 0.06; 0.06 0.15; 0.2 0.3];

figure;
sp_counter = 1;

for pen = 1:size(PN_units,1)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' PN_units{pen,1} '/tmp02/Spikes_' PN_units{pen,1} '_' PN_units{pen,2} '_Good_Pitch.mat']);

    PNUnits = PN_units{pen,3};
    Flist = unique(F0);
    repeats = unique(Y(:,5));

    for uu = 1:size(PNUnits,1)

        unit = PNUnits(uu,1);

        unitSpikes = Y(Y(:,3)==unit,:); % get spikes for just this unit
        profile = zeros(length(stims),17);

        if PNUnits(uu,2)==1
            window = [0 0.06];
        elseif PNUnits(uu,2)==2
            window = [0.06 0.15];
        else
            if pen<9
                window = [.3 .4];
            else
                window = [.2 .3];
            end

            continue
        end
    
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
            profile(ss,:) = zscore(meanSpikes);

        end % ends stim loop

        if sp_counter>25
            figure;
            sp_counter=1;
        end

        ax = subplot(5,5,sp_counter);
        imagesc(profile)

        yticks(1:length(stims));
        yticklabels(stims)
        xticks(1:4:17)
        xticklabels(Flist(1:4:17))
        sp_counter = sp_counter + 1;

    end % ends unit loop
end