%% Find tuning curves with correlations beyond chance using bootstrapping
% DEPENDENCIES: bootstrap_corr.m, plot_tuning_by_cond.m (if plot_yn == 'y')
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, January 2023

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

plot_yn = 'n'; % indicate whether you'd like to see plots of each tuning curve pairing and their correlation (will pause on each plot)

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

unit_corrs = cell(length(Animals),3); % allocate space to save the units we find as harmonicity neurons
ucount = 0;

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    stims_to_plot = unique(stims);
    stp_idxs = zeros(length(stims_to_plot));
    % find the indexes of the stims we want to evaluate
    for ss = 1:length(stims_to_plot)
        stp_idxs(ss) = find(strcmp(unique(type),stims_to_plot{ss}));
    end

    window = [0 0.1]; % response window

    uCorrs = zeros(size(BFs,1),size(BFs,2),13); % initialize space to save all of the correlations

    for uu = 1:length(units) % for each unit
        unit = units(uu);

        % now evaluate the correlation

        % evaluate every possible pairing of stimulus types
        for ss = 1:length(stims_to_plot)

            % only evaluate correlation if the neuron was sensitive to the F0 of this stimulus (ss)
            if sensitivity(uu,stp_idxs(ss))==1

                s = 1;

                while s < ss
                     % only evaluate correlation if the neuron was sensitive to the F0 of this stimulus (s)
                    if sensitivity(uu,stp_idxs(s))==1

                        [r,n5,n95] = bootstrap_corr(Y,type,F0,{stims_to_plot{ss},stims_to_plot{s}},unit);
                       
                        if r>n95 || r<n5 % if the correlation was outside the 95% CI of the null distribution
                            uCorrs(uu,stp_idxs(ss),stp_idxs(s)) = r;
                            uCorrs(uu,stp_idxs(s),stp_idxs(ss)) = r;
                            ucount = ucount + 1;
                
                            if strcmp(plot_yn,'y')
                
                                stim_to_plot = {stims_to_plot{ss},stims_to_plot{s}};
                                BFs_to_plot = [0 0]; % in order to not display best frequencies on the plot
                                plot_tuning_by_cond(Y,type,F0,unit,stim_to_plot,BFs_to_plot,Animals{ap},Pens{ap});
                                sgtitle('')
                                title(sprintf('r=%.2f, n95=%.2f, %s, %s, %d', r,n95,Animals{ap},Pens{ap},unit))
                
                                pause
                
                            end   
                        end
                    end

                    s = s + 1;

                end
            end
        end
    end % ends the unit loop

    % save all the correlation for each penetration
    unit_corrs{ap,1} = Animals{ap};
    unit_corrs{ap,2} = Pens{ap};
    unit_corrs{ap,3} = uCorrs;
    
end