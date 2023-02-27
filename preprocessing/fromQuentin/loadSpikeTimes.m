function [Y] = loadSpikeTimes(dataPath,qualityOpt,trigger_times,fs)
% Y : [spike absolute time - spike relative times - unit - stimulus # - repeat # - sweep #]
% Is the channel still relevant ? Or do we just work with units ?

if ~exist('dataPath','var') || isempty(dataPath)
    dataPath = cd;
end
if ~exist('qualityOpt','var') || isempty(qualityOpt)
    qualityOpt = 'Good';
end

% Loading infos
% load(fullfile(dataPath,'spike_times.mat'));
spike_times = readNPY(fullfile(dataPath,'spike_times.npy'));
load(fullfile(dataPath,'gridInfo.mat'));
% load(fullfile(dataPath,'expInfo.mat'));

spike_units = readNPY(fullfile(dataPath,'spike_clusters.npy'));
fid = fopen(fullfile(dataPath,'cluster_group.tsv'),'r');
dataArray = textscan(fid,'%f%s%[^\n\r]' , 'Delimiter', '\t', 'HeaderLines' ,1, 'ReturnOnError', false);
units = dataArray{1};
qualia = dataArray{2};

% Filter the units according the result of the manual clustering
switch qualityOpt 
    case 'Good'
        idx = strcmpi(qualia,'good');
    case 'All'
        idx =true(length(units),1);
    case 'MUA'
         idx = strcmpi(qualia,'mua');
    otherwise
        error('Unknown quality type.');
end
unitsList = units(idx);
idxSpikes = ismember(double(spike_units),unitsList);
spike_times = spike_times(idxSpikes);
spike_units = spike_units(idxSpikes);

% Build Y matrix
Y = zeros(sum(idxSpikes),6);
Y(:,1) = spike_times / fs * 1000; % spike times in ms
Y(:,3) = spike_units;

% sweep times
nSweepLoaded = length(trigger_times);
for i = 1:nSweepLoaded-1
%     idx = (Y(:,1) > (i-1)*expInfo.sweepLength) & (Y(:,1) < (i)*expInfo.sweepLength);
    idx = (Y(:,1) > trigger_times(i)) & (Y(:,1) < trigger_times(i+1));
    Y(idx,6) = i; % Sweep number
    Y(idx,4) = grid.randomisedGrid(i); % Stimulus number
    Y(idx,5) = sum(grid.randomisedGrid(1:i)==grid.randomisedGrid(i)); % Stimulus repeat number
    Y(idx,2) = Y(idx,1) - trigger_times(i); % Relative spike time
end
% Spike times in second
Y(:,[1 2]) = Y(:,[1 2]) /1000;

end