%% Edit spike mat to add stim info.

clear;
close all;

Animal = 'Derekah';
Pen = 'P01';
Qualia = 'Good';

addpath(genpath('C:\Users\Quentin Gaucher\code\npy-matlab'))

% dataPath = 'D:\Work\Data\ephys\P05-pitch2018\P05-pitch2018_sorted';
% dataPath = 'D:\Work\Data\ephys\Dorry\P04_Dory\P04-pitch2018\CRA';
dataPath = 'D:\data\pitch ephys\data\Derekah\P01\P01-pitch70dB2020_g0';
% load('D:\Work\Data\ephys\P05-pitch2018\P05-pitch2018_sorted\gridInfo.mat');
load(fullfile(dataPath,'gridInfo.mat'));
synch_ch = get_synch(dataPath);  % Get sweep start & stop

% min_trig_length_s = 0.7;
min_trig_length_s = 0.5;
% min_interTrig_length_s = 0.05;
min_interTrig_length_s = 0.12;
fs = 30000; %Hz

[start_time_ms] = get_triggers(synch_ch,min_trig_length_s,min_interTrig_length_s,fs); % Get triggers

% Y = get_spike_times(dataPath,'Good');
[Y] = loadSpikeTimes(dataPath,Qualia,start_time_ms,fs);
% Y -> [spike absolute time - spike relative times - unit - stimulus # - repeat # - sweep #]
% Grid info
stimNames = grid.stimFiles;
r = regexp(stimNames,'_(\d*)k?Hz','tokens');
for ii = 1:length(r)
    if isempty(r{ii})
        F0{ii} = '0';
    else
        F0{ii} = r{ii}{1}{1};
    end
end
F0 = F0';
F0 = str2double(F0);

type = regexp(stimNames,'_(\w*)_\d*k?Hz|_(\w*)_','tokens');
type = cellfun(@(x)((x{1})),type)';
type = cellfun(@(x)([x]),type,'UniformOutput',false);

if strcmp(Qualia,'MUA')
    save(['Spikes_' Animal '_' Pen '_' Qualia '_Pitch.mat'],'Y','type','F0','dataPath');
else
    save(['Spikes_' Animal '_' Pen '_Pitch.mat'],'Y','type','F0','dataPath');
end

disp('Done !')