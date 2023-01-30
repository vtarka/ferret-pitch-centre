
%% Script for evaluating the MI about pitch in the spike data as in Kerry's 2011 multiplexing paper
% DEPENDENCIES: MI2unbiased.m
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, January 2023


Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

rho_threshold = -1;
plot_yn = 'n'; % y = include plots of the unit's tuning, n = skip the plots

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

totalHN_count = 0;
unit_MIs  = cell(length(Animals),3); % allocate space to save the units we find as harmonicity neurons
rhos = [];

bin_dur = 30; % bin duration in ms
time_bins = 0:bin_dur/2:500;

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    allUnit_MI_tcs = zeros(length(units),length(stims),length(time_bins)-2);

    for uu = 1:length(units) % for each unit
        unit = units(uu);
        unitSpikes = Y(Y(:,3)==unit,:); % get the spikes of just this unit

        % pre-allocate space for the MI timecourses
        % for each stim we will generate a timecourse of MIs
        MI_timecourses = zeros(length(stims),length(time_bins)-2);

        % for each stimulus
        for ss = 1:length(stims)

            % for now only look at stims that modulated the neuronal
            % responses
            if sensitivity(uu,ss) ~= 1
                continue
            end

            for bin=1:length(time_bins)-2
                bin_start = time_bins(bin)/1000;
                bin_end = time_bins(bin+2)/1000;

                Spikes = zeros(length(repeats),length(Flist)); % create an emtpy matrix to hold spiking in each trial
               
                % for each F0
                for ff = 1:length(Flist)
            
                    stimNum = find(strcmp(type,stims{ss}) & (F0==Flist(ff)));  % unique name for combination of stim type and F0
        
                    if isempty(stimNum) % if this stim type and fo combo wasn't presented
                        nSpikes(:,ff) = 0;
                        continue
                    end
    
                    % for each presentation of this stim
                    for rr = 1:length(repeats)
                    
                        spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>=bin_start & unitSpikes(:,2)<bin_end;
                        nSpikes(rr,ff) = sum(spikeIDXs);
        
                    end   
                end % ends F0 loop

                uniqueSpikeRates = unique(nSpikes);

                % if the unit spiked the same rate all the time, the MI is
                % 0
                if length(uniqueSpikeRates)==1

                    MI_timecourses(ss,bin) = 0;

                else
                    
                    freq_histc = zeros(length(Flist),length(uniqueSpikeRates));
    
                    for ff = 1:length(Flist)
                        freq_histc(ff,:) = histcounts(nSpikes(:,ff),[uniqueSpikeRates; max(uniqueSpikeRates)+1]);
                    end
    
                    [~,rMI,bMI] = MI2unbiased(freq_histc);
                    tMI = rMI - bMI;
                    [aa,bb]=max(tMI(1:end-1)); %find the MI and bias we used
                    MI_timecourses(ss,bin) = rMI(bb);
                end

            end % ends bin loop
            
        end % ends stimulus loop

        allUnit_MI_tcs(uu,:,:) = MI_timecourses;

    end % ends unit loop

    % save all of this
    unit_MIs{ap,1} = Animals{ap};
    unit_MIs{ap,2} = Pens{ap};
    unit_MIs{ap,3} = allUnit_MI_tcs;

end % ends recording loop