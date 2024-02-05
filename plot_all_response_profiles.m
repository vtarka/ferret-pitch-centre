%% Plot all the response profiles in a single population
% DEPENDENCIES: get_response_profile.m
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, April 2023


mat_struct = load('HN_units_f.mat');
mat_cell = struct2cell(mat_struct);
units_by_rec = mat_cell{1};

stims = {'CT0','CT5','allHarm','low'};

window = [0 0.1]; % in seconds

figure;
uCounter = 1;

for pen = 1:length(units_by_rec)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/tmp02/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

    units = units_by_rec{pen,3};
%     units = unique(Y(:,3));

%     if pen < 9 
%         window = [0.3 0.4];
%     else
%         window = [0.2 0.3];
%     end

    for uu = 1:length(units)

%         if responsive(uu)==0 || sum(sensitivity(uu,:))==0
%             continue
%         end

        unit = units(uu);

        nexttile;
%         clf;    
        imagesc(get_response_profile(Y,type,F0,unit,stims,window))
        xticks([])
        yticks([])
        title(uCounter)
        uCounter = uCounter + 1;
%         yticks(1:length(stims))
%         yticklabels(stims)

%         pause

    end
end