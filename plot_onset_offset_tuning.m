% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

figure;

on_color = [0 0 0];
off_color = [1 0 0];

on_window = [0 0.1];

for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp02/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type); 
    if ap < 9
        off_window = [0.3 0.4];
    else
        off_window = [0.2 0.3];
    end
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    for uu = 1:length(units)
        if isempty(find(responsive_by_stim(uu,:), 1))
            continue
        end

        responsive_stims = find(responsive_by_stim(uu,:));
        for rs = 1:length(responsive_stims)
            clf;

            plot_tuning_by_cond(Y,type,F0,units(uu),stims(responsive_stims(rs)),0,Animals{ap},Pens{ap},on_window,on_color)
            plot_tuning_by_cond(Y,type,F0,units(uu),stims(responsive_stims(rs)),0,Animals{ap},Pens{ap},off_window,off_color)

            sgtitle('')
            sgtitle(stims(responsive_stims(rs)))

            pause

        end
    end
end