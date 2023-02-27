

%% Find HNs by defining them as neurons F0-sensitive to low harm and CT0 but not high harm
% low and CT0 tuning must also be correlated with rho > rho_threshold

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

totalonsetHN_count = 0;
totaloffsetHN_count = 0;
HN_units = cell(20,4); % allocate space to save the units we find as harmonicity neurons

timebins = 0:15:500;

uCounter = 1;

load('all_unit_MIs_30ms_2.mat')
load('all_peaks.mat')

% for each recording
for pen = 1:length(unit_MIs)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit_MIs{pen,1} '/tmp/Spikes_' unit_MIs{pen,1} '_' unit_MIs{pen,2} '_Good_Pitch.mat']);

    units = unique(Y(:,3));
    % want to find units where active(uu,1) & active(uu,11)

    onset_HNs = [];
    offset_HNs = [];
    for uu = 1:length(units)

        if all_peaks(uCounter,1,1) ~=0 && all_peaks(uCounter,11,1) ~= 0 && all_peaks(uCounter,10,1) == 0 % see if it has an onset response
            onset_HNs = [onset_HNs; units(uu)];
        end

        if all_peaks(uCounter,1,2) ~=0 && all_peaks(uCounter,11,2) ~= 0 && all_peaks(uCounter,10,2) == 0
            offset_HNs = [offset_HNs; units(uu)];
        end

        uCounter = uCounter + 1;

    end

    % save all the units we found for this penetration
    HN_units{pen,1} = unit_MIs{pen,1};
    HN_units{pen,2} = unit_MIs{pen,2};
    HN_units{pen,3} = onset_HNs;
    HN_units{pen,4} = offset_HNs;

end % ends recording loop


%% Plot the units we found

for pen = 1:length(unit_MIs)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit_MIs{pen,1} '/tmp/Spikes_' unit_MIs{pen,1} '_' unit_MIs{pen,2} '_Good_Pitch.mat']);

    onset_HNs = HN_units{pen,3};
    offset_HNs = HN_units{pen,4};

    stims = {'CT0','low','high'};

    for onHN = 1:length(onset_HNs)
        u = onset_HNs(onHN);
        plot_tuning_by_cond(Y,type,F0,u,stims,[0 0 0],HN_units{pen,1},HN_units{pen,2},[0 0.15])
    end

    for offHN = 1:length(offset_HNs)
        u = offset_HNs(offHN);
        if strcmp(HN_units{pen,1},'Noah')
            plot_tuning_by_cond(Y,type,F0,u,stims,[0 0 0],HN_units{pen,1},HN_units{pen,2},[0.3 0.45])
        else
            plot_tuning_by_cond(Y,type,F0,u,stims,[0 0 0],HN_units{pen,1},HN_units{pen,2},[0.2 0.35])
        end
    end

    pause
    close all

end % ends recording loop



%% plot histogram of distribution of MI values falling after onset vs offset
