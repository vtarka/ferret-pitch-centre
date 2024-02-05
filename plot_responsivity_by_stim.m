% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, May 2023

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

sound_responsive = 0;
not_sound_responsive = 0;
alpha = 0.05;

responsive_units_per_pen = zeros(length(Animals),2);

responsive_sums = zeros(13,1);

for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp02/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type); 
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    if ap < 9
        pre_window = [0.47 0.57];
        post_window = [0.3 0.4];
        responsive_by_stim(:,8:9) = [];
    else
        pre_window = [0.65 0.75];
        post_window = [0.2 0.3];
    end

    for ss = 1:13
        responsive_sums(ss) = responsive_sums(ss) + sum(responsive_by_stim(:,ss));
    end

end % ends loop through recordings


%% plot responsiveness by stimulus

figure;
bar(responsive_sums)
