
%% Get correlation between tuning curves of two conditions by taking 
% all possible pairings of trials 
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, January 2023

function [perm_corr, perm_std] = permutation_corr(Y,type,F0,stims,unit)

if length(stims) ~= 2
    return
end

repeats = unique(Y(:,5));
nRepeats = length(repeats);
Flist = unique(F0);
window = [0 0.1];

s1_s2_trials = zeros(2,nRepeats,length(Flist));
%s2_trials = zeros(nRepeats,length(Flist));

unitSpikes = Y(Y(:,3)==unit,:); % get the spikes of just this unit

% for each stim
for ss = 1:length(stims)
    
    nSpikes = zeros(length(repeats),length(Flist)); % create an emtpy matrix to hold spiking in each trial
    
    % for each F0
    for ff = 1:length(Flist)

        stimNum = find(strcmp(type,stims{ss}) & (F0==Flist(ff)));  % unique name for combination of stim type and F0

        if isempty(stimNum) % if this stim type and fo combo wasn't presented
            nSpikes(:,ff) = 0;
            continue
        end

        % for each presentation of this stim
        for rr = 1:length(repeats)
        
            spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
            nSpikes(rr,ff) = sum(spikeIDXs);

            s1_s2_trials(ss,rr,ff) = sum(spikeIDXs)/diff(window);
        end   
    end % ends F0 loop

end % end stim loop

s1_trials = squeeze(s1_s2_trials(1,:,:));
s2_trials = squeeze(s1_s2_trials(2,:,:));

corr_dist = zeros(size(s1_trials,1)*size(s2_trials,1),1);
corr_counter = 1;

for s1 = 1:size(s1_trials,1)
    for s2 = 1:size(s2_trials,1)
        corr_dist(corr_counter) = corr(s1_trials(s1,:)',s2_trials(s2,:)');
        corr_counter = corr_counter + 1;
    end
end

perm_corr = mean(corr_dist,'omitnan');
perm_std = std(corr_dist,'omitnan');

end