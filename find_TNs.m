%% Find temporal neurons by F0-sensitivity to high, low, and CT0
% DEPENDENCIES: plot_tuning_by_cond.m (if plot_yn == 'y')
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, March 2023

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

locs = [1 1 2 1 5 3 5 4 3 3 1 2 3 2 1 1 1 3 4 3]; % a list of locations by penetration (aligns with Animals and Pens variables)
loc_counts = [];

% % location list:    low A1      high A1    low AAF     high AAF    PPF
% %                      1          2           3            4        5

plot_yn = 'n'; % y = include plots of the unit's tuning, n = skip the plots

shuffle_tuning_yn = 'y'; % y = shuffle the tuning profiles to see if the unit is more aligned than random chance, n = skip this
null_percentile_threshold = 99; % specify which percentile of the null distribution to use as the threshold for significant tuning alignment (only used if shuffle_tuning_yn == 'y')
nNullRuns = 1000; % specify how many times to shuffle the tuning to get the null distribution (only used if shuffle_tuning_yn == 'y')

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

totalTN_count = 0;
TN_units = cell(length(Animals),3); % allocate space to save the units we find as harmonicity neurons

stims_for_profile = {'CT0','CT5','CT10','allHarm','high'};
stims_to_plot = {'CT0','low','high'};

window = [0 0.1];

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp02/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    units = unique(Y(:,3));

    TN_unit_list = []; % to keep track of TNs we find
    for uu = 1:length(units) % for each unit
        unit = units(uu);

        high_stim_num = find(strcmp(stims,'high'));
        CT0_stim_num = 1; % it's always 1

        % if the unit is F0-sensitive to CT0 and high
        if sensitivity(uu,CT0_stim_num) && sensitivity(uu,high_stim_num)

            if strcmp(shuffle_tuning_yn,'y')

                profile = get_response_profile(Y,type,F0,unit,stims_for_profile,window);

                real_corr = get_avg_pairwise_corr(profile);

                shuffled_profile = profile;
                shuffled_corrs = zeros(nNullRuns,1);
                for i = 1:nNullRuns
                    for ss = 1:length(stims_for_profile)
                        
                        shuffle_key = randperm(size(profile,2)); % random permutation of the number of frequencies we presented

                        shuffled_profile(ss,:) = profile(ss,shuffle_key);
                    end

                    shuffled_corrs(i) = get_avg_pairwise_corr(shuffled_profile);
                end

                null_grand_avg_corr = prctile(shuffled_corrs,null_percentile_threshold);

                % if our real correlation (real tuning alignment) exceeds the threshold in the null distribution
                if real_corr > null_grand_avg_corr
                    loc_counts(end+1) = locs(ap); % we've found an TN
                    totalTN_count = totalTN_count + 1;
                    TN_unit_list = [TN_unit_list; unit];

                    if strcmp(plot_yn,'y')

                        clf
                        BFs_to_plot = [0 0 0]; %BFs(uu,[10 11 1]);
                        plot_tuning_by_cond(Y,type,F0,units(uu),stims_to_plot,BFs_to_plot,Animals{ap},Pens{ap},window);
         
                        pause
                    end  
                end
            else

                loc_counts(end+1) = locs(ap); % we've found an TN
                totalTN_count = totalTN_count + 1;
                TN_unit_list = [TN_unit_list; unit];

                if strcmp(plot_yn,'y')

                    clf
                    BFs_to_plot = [0 0 0];
                    plot_tuning_by_cond(Y,type,F0,units(uu),stims_to_plot,BFs_to_plot,Animals{ap},Pens{ap},window);
     
                    pause
                end  
            end % ends shuffle tuning if case
        end % ends F0-sensitivity if case
    end % ends the unit loop

    % save all the units we found for this penetration
    TN_units{ap,1} = Animals{ap};
    TN_units{ap,2} = Pens{ap};
    TN_units{ap,3} = TN_unit_list;

end % ends recording loop