
Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

locs = [1 1 2 1 5 3 5 4 4 1 1 1 1 2 1 1 1 4 2 3];

sensitive_by_region = zeros(5,1);
units_by_region = zeros(5,1);

for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/Batch0225/Spikes_' Animals{ap} '_' Pens{ap} '_' 'Good' '_Pitch.mat']);

    units = unique(Y(:,3));

    nunits = length(units);

%     sensitive_units = find(sensitivity==1);
    nsensitive = length(find(sensitivity(:,1)==1));

    sensitive_by_region(locs(ap)) = sensitive_by_region(locs(ap)) + nsensitive;
    units_by_region(locs(ap)) = units_by_region(locs(ap)) + nunits;

    sprintf('%s %s %d total and %d sensitive',Animals{ap},Pens{ap},nunits,nsensitive)
end