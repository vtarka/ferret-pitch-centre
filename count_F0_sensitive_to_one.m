

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

window = [0 0.1]; % in ms
not_sensitive_count = 0;

% for each recording
for ap = 1:length(Animals)


    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp02/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);
    units = unique(Y(:,3));

    for uu = 1:length(units)

        if responsive(uu)

            if isempty(find(sensitivity(uu,:)))
                not_sensitive_count = not_sensitive_count + 1;
            end
        end
    end

end