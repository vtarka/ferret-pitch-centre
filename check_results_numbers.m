Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

region = [1 1 2 1 5 3 5 4 3 3 1 2 3 2 1 1 1 3 4 3];
n_neurons_by_region = zeros(length(unique(region)),1);
n_neurons = 0;
n_responsive_neurons = 0;
n_sensitive_neurons = 0;

sensitivity_counts = [];

sensitive_stims = [];

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/Batch0824/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type); 
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    n_neurons = n_neurons + length(units);
    n_neurons_by_region(region(ap)) = n_neurons_by_region(region(ap)) + length(units);

    n_responsive_neurons = n_responsive_neurons + length(find(responsive));

    if ap<9
        pre_sensitivity(:,[8 10])= [];
    end
    
    summed_sensitivity = sum(pre_sensitivity,2);
    n_sensitive_neurons = n_sensitive_neurons + length(find(summed_sensitivity));

    sensitivity_counts = [sensitivity_counts; summed_sensitivity(summed_sensitivity~=0)];


    [~,c] = find(pre_sensitivity);
    sensitive_stims = [sensitive_stims; c];

end