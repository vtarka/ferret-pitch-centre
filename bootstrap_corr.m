
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
window = [0 0.1];

s1_s2_trials = zeros(2,nRepeats,length(Flist));

unitSpikes = Y(Y(:,3)==unit,:); % get the spikes of just this unit

tuning = zeros(2,17);

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

            s1_s2_trials(ss,rr,ff) = sum(spikeIDXs);%/diff(window);
        end   
    end % ends F0 loop

    tuning(ss,:) = mean(s1_s2_trials(ss,:,:));

end % end stim loop

s1_trials = squeeze(s1_s2_trials(1,:,:));
s2_trials = squeeze(s1_s2_trials(2,:,:));

if strcmp('low',stims{1}) || strcmp('low',stims{2})
    s1_trials(:,1:2) = [];
    s2_trials(:,1:2) = [];
    Flist(1:2) = [];
end

s1_rep = [];
s2_rep = [];

for f0 = 1:length(Flist)

    idx = randperm(13*13);
    s1 = reshape(repmat(s1_trials(:,f0)',13,1),[],1);
    s2 = repmat(s2_trials(:,f0),13,1);    
    
    s1_rep = [s1_rep; s1];
    s2_rep = [s2_rep; s2];

end

r = corr(s1_rep,s2_rep);

% figure('Position',[1900 500 1800 1200])
% subplot(3,1,1); scatter(1:length(s1_rep),s1_rep,70,'MarkerFaceColor','k','MarkerFaceAlpha',0.2)
% hold on; scatter(1:length(s2_rep),s2_rep,70,'MarkerFaceColor','r','MarkerFaceAlpha',0.2)
% xticks(0:13*13:13*13*17)
% title(r)
% 
% % subplot(2,1,2); scatter(1:length(s2_rep),s2_rep,'filled')
% % xticks(0:13*13:13*13*17)
% 
% nonrep_s1 = reshape(s1_trials,[],1);
% nonrep_s2 = reshape(s2_trials,[],1);
% nrr = corr(nonrep_s1,nonrep_s2);
% subplot(3,1,2)
% scatter(1:length(nonrep_s1),nonrep_s1,70,'MarkerFaceColor','k','MarkerFaceAlpha',0.2)
% hold on; scatter(1:length(nonrep_s2),nonrep_s2,70,'MarkerFaceColor','r','MarkerFaceAlpha',0.2)
% title(nrr)
% 
% tr = corr(tuning(1,:)',tuning(2,:)');
% subplot(3,1,3)
% plot(tuning(1,:),'k')
% hold on; plot(tuning(2,:),'r')
% title(tr)
% 
% pause

null_perms = 1000;
null_rhos = zeros(null_perms,1);

for p = 1:null_perms

% %     shuffled_s1 = s1_trials(randperm(m*n));
    shuffled_s2_rep = s2_rep(randperm(length(s2_rep)));
% 
% %     ss1_vector = reshape(shuffled_s1,[],1);
% %     ss2_vector = reshape(shuffled_s2,[],1);

    null_rhos(p) = corr(s1_rep,shuffled_s2_rep);
end

null_5 = prctile(null_rhos,5);
null_95 = prctile(null_rhos,95);

end