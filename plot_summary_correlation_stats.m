
%% Find average correlations 
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, January 2023

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

all_corrs = zeros(839,13,13);
i=1;

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    window = [0 0.15]; % response window

    HN_unit_list = []; % to keep track of HNs we find
    for uu = 1:length(units) % for each unit
       unit = units(uu);

       % if this unit isn't F0-sensitive to any stimulus 
       if isempty(find(sensitivity(uu,:), 1))
           continue
       end

       unitSpikes = Y(Y(:,3)==units(uu),:); % spikes for just this unit

       tuning = zeros(length(stims),length(Flist));

        % for each stim type
       for ss = 1:length(stims)

            spike_counts = zeros(length(repeats),length(Flist)); % initialize space to save the number of spikes evoked
    
            % for each F0
            for ff = 1:length(Flist)
    
                stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); % unique name for combination of stim type and F0

                if isempty(stimNum) % if this stim type and F0 combo wasn't presented
                    continue
                end
    
                % for each presentation of this stim
                for rr = 1:length(repeats)
                
                    spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                    nSpikes = sum(spikeIDXs);
    
                    spike_counts(rr,ff) = nSpikes ./ diff(window); % convert to spike rate and save
    
                end 
            end

            meanSpikes = mean(spike_counts);
            tuning(ss,:) = meanSpikes;
       end

        rho_matrix = corrcoef(tuning');
    
        all_corrs(i,:,:) = rho_matrix;
        i = i+1;
    end
end

mean_corrs = squeeze(mean(all_corrs,'omitnan'));
mean_corrs = mean_corrs - diag(diag(mean_corrs));
figure; imagesc(mean_corrs); colorbar
xticks(1:13); xticklabels(stims)
yticks(1:13); yticklabels(stims)


%% Find average correlations ignoring stimuli that aren't frequency sensitive

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

all_corrs = zeros(839,13,13);
i=1;

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    window = [0 0.15]; % response window

    HN_unit_list = []; % to keep track of HNs we find
    for uu = 1:length(units) % for each unit
       unit = units(uu);

       % if this unit isn't F0-sensitive to any stimulus 
       if isempty(find(sensitivity(uu,:), 1))
           continue
       end

       unitSpikes = Y(Y(:,3)==units(uu),:); % spikes for just this unit

       tuning = zeros(length(stims),length(Flist));

        % for each stim type
       for ss = 1:length(stims)
    
            spike_counts = zeros(length(repeats),length(Flist)); % initialize space to save the number of spikes evoked
    
            % for each F0
            for ff = 1:length(Flist)
    
                stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); % unique name for combination of stim type and F0

                if isempty(stimNum) % if this stim type and F0 combo wasn't presented
                    continue
                end
    
                % for each presentation of this stim
                for rr = 1:length(repeats)
                
                    spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                    nSpikes = sum(spikeIDXs);
    
                    spike_counts(rr,ff) = nSpikes ./ diff(window); % convert to spike rate and save
    
                end 
            end

            meanSpikes = mean(spike_counts);
            tuning(ss,:) = meanSpikes;
       end

        rho_matrix = corrcoef(tuning');
        not_F0_sensitive = find(sensitivity(uu,:)==0);
        rho_matrix(not_F0_sensitive,:) = nan;
        rho_matrix(:,not_F0_sensitive,:) = nan;
    
        all_corrs(i,:,:) = rho_matrix;
        i = i+1;
    end
end

mean_corrs = squeeze(mean(all_corrs,'omitnan'));
mean_corrs = mean_corrs - diag(diag(mean_corrs));
figure; imagesc(mean_corrs); a = colorbar;
xticks(1:13); xticklabels(stims)
yticks(1:13); yticklabels(stims)
set(gca,'fontsize',22)
ylabel(a,'Avg correlation','Fontsize',26)
title('Avg correlations of F0-sensitive tuning')


%% Evaluate average correlation by each anatomical region

%%%%%%%%%%%%%%%%%%%%% high/mid A1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Animals = {'Noah','Ronnie','Derry'};

Pens = {'P03','P13','P03'};

Qualia = 'Good';

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

all_corrs = zeros(839,13,13);
i=1;

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    window = [0 0.15]; % response window

    HN_unit_list = []; % to keep track of HNs we find
    for uu = 1:length(units) % for each unit
       unit = units(uu);

       % if this unit isn't F0-sensitive to any stimulus 
       if isempty(find(sensitivity(uu,:), 1))
           continue
       end

       unitSpikes = Y(Y(:,3)==units(uu),:); % spikes for just this unit

       tuning = zeros(length(stims),length(Flist));

        % for each stim type
       for ss = 1:length(stims)
    
            spike_counts = zeros(length(repeats),length(Flist)); % initialize space to save the number of spikes evoked
    
            % for each F0
            for ff = 1:length(Flist)
    
                stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); % unique name for combination of stim type and F0

                if isempty(stimNum) % if this stim type and F0 combo wasn't presented
                    continue
                end
    
                % for each presentation of this stim
                for rr = 1:length(repeats)
                
                    spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                    nSpikes = sum(spikeIDXs);
    
                    spike_counts(rr,ff) = nSpikes ./ diff(window); % convert to spike rate and save
    
                end 
            end

            meanSpikes = mean(spike_counts);
            tuning(ss,:) = meanSpikes;
       end

        rho_matrix = corrcoef(tuning');
        not_F0_sensitive = find(sensitivity(uu,:)==0);
        rho_matrix(not_F0_sensitive,:) = nan;
        rho_matrix(:,not_F0_sensitive,:) = nan;
    
        all_corrs(i,:,:) = rho_matrix;
        i = i+1;
    end
end

all_corrs(i:end,:,:) = [];

mean_corrs = squeeze(mean(all_corrs,'omitnan'));
mean_corrs = mean_corrs - diag(diag(mean_corrs));
figure; imagesc(mean_corrs); a = colorbar;
xticks(1:13); xticklabels(stims)
yticks(1:13); yticklabels(stims)
set(gca,'fontsize',22)
ylabel(a,'Avg correlation','Fontsize',26)
title('Avg correlations of F0-sensitive tuning in high/mid A1')

%%
%%%%%%%%%%%%%%%%%%%%% low A1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Animals = {'Noah','Noah','Noah','Ronnie','Derry','Dory'};

Pens = {'P01','P02','P04','P08','P05','P00'};

Qualia = 'Good';

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

all_corrs = zeros(839,13,13);
i=1;

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    window = [0 0.15]; % response window

    HN_unit_list = []; % to keep track of HNs we find
    for uu = 1:length(units) % for each unit
       unit = units(uu);

       % if this unit isn't F0-sensitive to any stimulus 
       if isempty(find(sensitivity(uu,:), 1))
           continue
       end

       unitSpikes = Y(Y(:,3)==units(uu),:); % spikes for just this unit

       tuning = zeros(length(stims),length(Flist));

        % for each stim type
       for ss = 1:length(stims)
    
            spike_counts = zeros(length(repeats),length(Flist)); % initialize space to save the number of spikes evoked
    
            % for each F0
            for ff = 1:length(Flist)
    
                stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); % unique name for combination of stim type and F0

                if isempty(stimNum) % if this stim type and F0 combo wasn't presented
                    continue
                end
    
                % for each presentation of this stim
                for rr = 1:length(repeats)
                
                    spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                    nSpikes = sum(spikeIDXs);
    
                    spike_counts(rr,ff) = nSpikes ./ diff(window); % convert to spike rate and save
    
                end 
            end

            meanSpikes = mean(spike_counts);
            tuning(ss,:) = meanSpikes;
       end

        rho_matrix = corrcoef(tuning');
        not_F0_sensitive = find(sensitivity(uu,:)==0);
        rho_matrix(not_F0_sensitive,:) = nan;
        rho_matrix(:,not_F0_sensitive,:) = nan;
    
        all_corrs(i,:,:) = rho_matrix;
        i = i+1;
    end
end

all_corrs(i:end,:,:) = [];

mean_corrs = squeeze(mean(all_corrs,'omitnan'));
mean_corrs = mean_corrs - diag(diag(mean_corrs));
figure; imagesc(mean_corrs); a = colorbar;
xticks(1:13); xticklabels(stims)
yticks(1:13); yticklabels(stims)
set(gca,'fontsize',22)
ylabel(a,'Avg correlation','Fontsize',26)
title('Avg correlations of F0-sensitive tuning in low A1')

%%
%%%%%%%%%%%%%%%%%%%%% low AAF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Animals = {'Noah','Ronnie','Ronnie','Derry','Dory','Dory'};

Pens = {'P06','P04','P05','P02','P01','P04'};

Qualia = 'Good';

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

all_corrs = zeros(839,13,13);
i=1;

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    window = [0 0.15]; % response window

    HN_unit_list = []; % to keep track of HNs we find
    for uu = 1:length(units) % for each unit
       unit = units(uu);

       % if this unit isn't F0-sensitive to any stimulus 
       if isempty(find(sensitivity(uu,:), 1))
           continue
       end

       unitSpikes = Y(Y(:,3)==units(uu),:); % spikes for just this unit

       tuning = zeros(length(stims),length(Flist));

        % for each stim type
       for ss = 1:length(stims)
    
            spike_counts = zeros(length(repeats),length(Flist)); % initialize space to save the number of spikes evoked
    
            % for each F0
            for ff = 1:length(Flist)
    
                stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); % unique name for combination of stim type and F0

                if isempty(stimNum) % if this stim type and F0 combo wasn't presented
                    continue
                end
    
                % for each presentation of this stim
                for rr = 1:length(repeats)
                
                    spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                    nSpikes = sum(spikeIDXs);
    
                    spike_counts(rr,ff) = nSpikes ./ diff(window); % convert to spike rate and save
    
                end 
            end

            meanSpikes = mean(spike_counts);
            tuning(ss,:) = meanSpikes;
       end

        rho_matrix = corrcoef(tuning');
        not_F0_sensitive = find(sensitivity(uu,:)==0);
        rho_matrix(not_F0_sensitive,:) = nan;
        rho_matrix(:,not_F0_sensitive,:) = nan;
    
        all_corrs(i,:,:) = rho_matrix;
        i = i+1;
    end
end

all_corrs(i:end,:,:) = [];

mean_corrs = squeeze(mean(all_corrs,'omitnan'));
mean_corrs = mean_corrs - diag(diag(mean_corrs));
figure; imagesc(mean_corrs); a = colorbar;
xticks(1:13); xticklabels(stims)
yticks(1:13); yticklabels(stims)
set(gca,'fontsize',22)
ylabel(a,'Avg correlation','Fontsize',26)
title('Avg correlations of F0-sensitive tuning in low AAF')


%%
%%%%%%%%%%%%%%%%%%%%% high AAF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Animals = {'Noah'};

Pens = {'P08'};

Qualia = 'Good';

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

all_corrs = zeros(839,13,13);
i=1;

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    window = [0 0.15]; % response window

    HN_unit_list = []; % to keep track of HNs we find
    for uu = 1:length(units) % for each unit
       unit = units(uu);

       % if this unit isn't F0-sensitive to any stimulus 
       if isempty(find(sensitivity(uu,:), 1))
           continue
       end

       unitSpikes = Y(Y(:,3)==units(uu),:); % spikes for just this unit

       tuning = zeros(length(stims),length(Flist));

        % for each stim type
       for ss = 1:length(stims)
    
            spike_counts = zeros(length(repeats),length(Flist)); % initialize space to save the number of spikes evoked
    
            % for each F0
            for ff = 1:length(Flist)
    
                stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); % unique name for combination of stim type and F0

                if isempty(stimNum) % if this stim type and F0 combo wasn't presented
                    continue
                end
    
                % for each presentation of this stim
                for rr = 1:length(repeats)
                
                    spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                    nSpikes = sum(spikeIDXs);
    
                    spike_counts(rr,ff) = nSpikes ./ diff(window); % convert to spike rate and save
    
                end 
            end

            meanSpikes = mean(spike_counts);
            tuning(ss,:) = meanSpikes;
       end

        rho_matrix = corrcoef(tuning');
        not_F0_sensitive = find(sensitivity(uu,:)==0);
        rho_matrix(not_F0_sensitive,:) = nan;
        rho_matrix(:,not_F0_sensitive,:) = nan;
    
        all_corrs(i,:,:) = rho_matrix;
        i = i+1;
    end
end

all_corrs(i:end,:,:) = [];

mean_corrs = squeeze(mean(all_corrs,'omitnan'));
mean_corrs = mean_corrs - diag(diag(mean_corrs));
figure; imagesc(mean_corrs); a = colorbar;
xticks(1:13); xticklabels(stims)
yticks(1:13); yticklabels(stims)
set(gca,'fontsize',22)
ylabel(a,'Avg correlation','Fontsize',26)
title('Avg correlations of F0-sensitive tuning in high AAF')

%%
%%%%%%%%%%%%%%%%%%%%% low PPF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Animals = {'Noah','Noah'};

Pens = {'P05','P07'};

Qualia = 'Good';

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

all_corrs = zeros(839,13,13);
i=1;

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    window = [0 0.15]; % response window

    HN_unit_list = []; % to keep track of HNs we find
    for uu = 1:length(units) % for each unit
       unit = units(uu);

       % if this unit isn't F0-sensitive to any stimulus 
       if isempty(find(sensitivity(uu,:), 1))
           continue
       end

       unitSpikes = Y(Y(:,3)==units(uu),:); % spikes for just this unit

       tuning = zeros(length(stims),length(Flist));

        % for each stim type
       for ss = 1:length(stims)
    
            spike_counts = zeros(length(repeats),length(Flist)); % initialize space to save the number of spikes evoked
    
            % for each F0
            for ff = 1:length(Flist)
    
                stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); % unique name for combination of stim type and F0

                if isempty(stimNum) % if this stim type and F0 combo wasn't presented
                    continue
                end
    
                % for each presentation of this stim
                for rr = 1:length(repeats)
                
                    spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                    nSpikes = sum(spikeIDXs);
    
                    spike_counts(rr,ff) = nSpikes ./ diff(window); % convert to spike rate and save
    
                end 
            end

            meanSpikes = mean(spike_counts);
            tuning(ss,:) = meanSpikes;
       end

        rho_matrix = corrcoef(tuning');
        not_F0_sensitive = find(sensitivity(uu,:)==0);
        rho_matrix(not_F0_sensitive,:) = nan;
        rho_matrix(:,not_F0_sensitive,:) = nan;
    
        all_corrs(i,:,:) = rho_matrix;
        i = i+1;
    end
end

all_corrs(i:end,:,:) = [];

mean_corrs = squeeze(mean(all_corrs,'omitnan'));
mean_corrs = mean_corrs - diag(diag(mean_corrs));
figure; imagesc(mean_corrs); a = colorbar;
xticks(1:13); xticklabels(stims)
yticks(1:13); yticklabels(stims)
set(gca,'fontsize',22)
ylabel(a,'Avg correlation','Fontsize',26)
title('Avg correlations of F0-sensitive tuning in low PPF')