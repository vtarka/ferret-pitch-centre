
%% Determine whether a neuron was active/responsive or not based on its MI timecourse
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, February 2023

% max_MIs = [];
for pen = 1:length(unit_MIs)

    MI_tcs = unit_MIs{pen,3};

%     max_MIs = [max_MIs; squeeze(max(MI_tcs,[],3))];

    max_MIs = squeeze(max(MI_tcs,[],3));
    active = max_MIs > 0.12;

    unit_MIs{pen,4} = active;

end % ends recording loop

save('all_unit_MIs_30ms_2.mat','unit_MIs')

% figure;
% histogram(max_MIs)


% threshold at 0.07 for 20 ms