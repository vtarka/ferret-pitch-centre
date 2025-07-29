%% Find pitch neurons by F0-sensitivity to high, low, and CT0
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
null_percentile_threshold = 95; % specify which percentile of the null distribution to use as the threshold for significant tuning alignment (only used if shuffle_tuning_yn == 'y')
nNullRuns = 10000; % specify how many times to shuffle the tuning to get the null distribution (only used if shuffle_tuning_yn == 'y')
nullDistributions = zeros(50,nNullRuns);
uCounter = 1;

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

totalPN_count = 0;
total_discarded = 0;
PN_units = cell(length(Animals),4); % allocate space to save the units we find as harmonicity neurons

stims_for_profile = {'CT0','CT5','allHarm','high','low'};
stims_to_plot = {'CT0','low','high'};

window = [0 0.1];

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/Batch0225/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    units = unique(Y(:,3));

    PN_unit_list = []; % to keep track of PNs we find
    discarded_unit_list = [];

    for uu = 1:length(units) % for each unit
        unit = units(uu);

        high_stim_num = find(strcmp(stims,'high'));
        low_stim_num = find(strcmp(stims,'low'));
        CT0_stim_num = find(strcmp(stims,'CT0')); % it's always 1

        % if the unit is F0-sensitive to CT0 and high
        if sensitivity(uu,high_stim_num) && sensitivity(uu,low_stim_num) && sensitivity(uu,CT0_stim_num)

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
                    loc_counts(end+1) = locs(ap); % we've found an PN
                    totalPN_count = totalPN_count + 1;
                    PN_unit_list = [PN_unit_list; unit];
                    nullDistributions(uCounter,:) = shuffled_corrs;
                    uCounter = uCounter + 1;

                    if strcmp(plot_yn,'y')

                        clf
                        BFs_to_plot = [0 0 0]; %BFs(uu,[10 11 1]);
                        plot_tuning_by_cond(Y,type,F0,unit,stims_to_plot,BFs_to_plot,Animals{ap},Pens{ap},window);
         
                        pause
                    end 

%                     figure(1);
%                     nexttile;
%                     imagesc(get_response_profile(Y,type,F0,unit,stims_for_profile,window))
%                     xticks([])
%                     yticks([])
                else
                    total_discarded = total_discarded + 1;
%                     figure(2);
%                     nexttile;
%                     imagesc(get_response_profile(Y,type,F0,unit,stims_for_profile,window))
%                     xticks([])
%                     yticks([])

                    discarded_unit_list = [discarded_unit_list; unit];
                end
            else

                loc_counts(end+1) = locs(ap); % we've found an PN
                totalPN_count = totalPN_count + 1;
                PN_unit_list = [PN_unit_list; unit];

                if strcmp(plot_yn,'y')

                    clf
                    BFs_to_plot = [0 0 0];
                    plot_tuning_by_cond(Y,type,F0,unit,stims_to_plot,BFs_to_plot,Animals{ap},Pens{ap},window);
     
                    pause
                end  
            end % ends shuffle tuning if case
        end % ends F0-sensitivity if case
    end % ends the unit loop

    % save all the units we found for this penetration
    PN_units{ap,1} = Animals{ap};
    PN_units{ap,2} = Pens{ap};
    PN_units{ap,3} = PN_unit_list;
    PN_units{ap,4} = discarded_unit_list;

    nullDistributions(totalPN_count+1:end,:) = [];

end % ends recording loop

figure(2); sgtitle('Discarded Units')
figure(1); sgtitle('Pitch Units')
% 
% save('PN_units','PN_units')
% save('PN_nulls','nullDistributions')