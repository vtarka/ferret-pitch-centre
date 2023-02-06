
%% Find neurons with a MI peak in response to stim onset and stim offset
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, February 2023

load('unit_MIs_30ms.mat')

twoPeak_units = cell(20,3); % allocate space to save the units we find as harmonicity neurons

for pen = 1:length(unit_MIs)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' uMI50ms{pen,1} '/tmp/Spikes_' uMI50ms{pen,1} '_' uMI50ms{pen,2} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    HN_unit_list = []; % to keep track of HNs we find
    for uu = 1:length(units) % for each unit
    end
    
end % ends recording  loop