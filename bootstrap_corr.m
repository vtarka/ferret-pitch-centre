
%% Get correlation between tuning curves of two conditions by taking 
% all possible pairings of trials 
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, January 2023

function [r,null_5,null_95] = bootstrap_corr(Y,type,F0,stims,unit)

if length(stims) ~= 2
    return
end

repeats = unique(Y(:,5));
nRepeats = length(repeats);
Flist = unique(F0);
window = [0 0.08];

s1_s2_trials = zeros(2,nRepeats,length(Flist));

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

s1_rep = [];
s2_rep = [];

for f0 = 1:length(Flist)

    idx = randperm(13*13);
    s1 = reshape(repmat(s1_trials(:,f0)',13,1),[],1);
    s2 = repmat(s2_trials(:,f0),13,1);    
    
    s1_rep = [s1_rep; s1];
    s2_rep = [s2_rep; s2];

end

figure; subplot(2,1,1); scatter(1:length(s1_rep),s1_rep,'filled')
xticks(0:13*13:13*13*17)

subplot(2,1,2); scatter(1:length(s2_rep),s2_rep,'filled')
xticks(0:13*13:13*13*17)

r = corr(s1_rep,s2_rep);

null_perms = 10000;
null_rhos = zeros(null_perms,1);

for p = 1:null_perms

%     shuffled_s1 = s1_trials(randperm(m*n));
    shuffled_s2_rep = s2_rep(randperm(length(s2_rep)));

%     ss1_vector = reshape(shuffled_s1,[],1);
%     ss2_vector = reshape(shuffled_s2,[],1);

    null_rhos(p) = corr(s1_rep,shuffled_s2_rep);
end

null_5 = prctile(null_rhos,5);
null_95 = prctile(null_rhos,95);

end