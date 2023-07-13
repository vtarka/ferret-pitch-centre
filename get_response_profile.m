%% Function to generate the response profile for a single neuron
% A response profile is an nStim by nPitch matrix where each row is the z-scored avg tuning curve for that stimulus
% DEPNDENCIES: none
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, April 2023


function profile = get_response_profile(Y,type,F0,unit,stims_for_profile,window)

% INPUTS:
% Y - n x 6 matrix where Y(:,)=spiketimes, Y(:,3)=unit spiking, Y(:,4)=stimulus, Y(:,5)=repeat, Y(:,6)=trial number
% type - nStimuli x 1 vector labeling the sound type (CT0, tone, etc) of each stimulus
% F0 - nStimuli x 1 vector labeling the F0 of each stimulus
% unit - the unit to be plotted
% stims_for_profile - cell array of strings containing stimuli names to be included in profile
% window - vector with two elements: window(1) is the start of the response window, window(2) is the end (in ms)


Flist = unique(F0);
repeats = unique(Y(:,5));

unitSpikes = Y(Y(:,3)==unit,:); % get the spikes of just this unit
profile = zeros(length(stims_for_profile),17); % allocate space to construct the tuning profile

% go through each stim we want to plot
for ss = 1:length(stims_for_profile)

    nSpikes = zeros(length(repeats),length(Flist)); % allocate space to save spiking info for each trial
        
    for ff = 1:length(Flist) % for each frequency

        stimNum = find(strcmp(type,stims_for_profile{ss}) & (F0==Flist(ff))); % unique name for combination of stim type and F0

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
    z_spikes = zscore(meanSpikes);
    norm_z_spikes = (z_spikes - min(z_spikes)) / (max(z_spikes) - min(z_spikes));
    profile(ss,:) = z_spikes; % save the z-scored tuning curve to the profile

    % the first two pitches were not presented for the low harmonics 
    if strcmp('low',stims_for_profile{ss}) || strcmp('F0MaskLow',stims_for_profile{ss})
        profile(ss,1:2) = nan;
    end

end % ends stim loop
end