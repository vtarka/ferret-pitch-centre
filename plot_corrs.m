
%% Quick plot to test out correlation assesment methods
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, January 2023


figure('Position',[1900 500 1800 1200])
stims = {'low','CT0'};
colors = colormap(hsv(length(stims))); % make the colormap to be used later
sp = 0;

% for each penetration
for pen = 1:length(HN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' HN_units{pen,1} '/tmp/Spikes_' HN_units{pen,1} '_' HN_units{pen,2} '_Good_Pitch.mat'])

    stims = {'low','CT0'};
    Flist = unique(F0);
%     repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

%     window = [0 0.15]; % response window

    HNUnits = HN_units{pen,3}; % array of units we found to be harmonicity neurons
    [~,HNUnit_IDXs] = ismember(HNUnits,allUnits); % find the indices of these pitch neurons within the list of all units

    HN_BFs = BFs(HNUnit_IDXs,:); % get the PN's best frequencies

    for hn = 1:length(HNUnits)
        
        [r, n5, n95] = bootstrap_corr(Y,type,F0,stims,HNUnits(hn));

        plot_tuning_by_cond(Y,type,F0,HNUnits(hn),stims,HN_BFs(hn,[11 1]),HN_units{pen,1},HN_units{pen,2});
        sgtitle('');
        title(sprintf('RHO: %.2f, [5 95]: [%.2f  %.2f]',r,n5,n95))
        pause
    end

end



%%

% load('all_unit_corrs.mat')

all_corrs = [];

for ap = 1:length(unit_corrs)
    all_corrs = [all_corrs; unit_corrs{ap,3}];
end

all_corrs_vec = reshape(all_corrs,[],1);
all_corrs_vec(all_corrs_vec==0) = [];

figure; histogram(all_corrs_vec)