%% Function to measure how sensitive a neuron is to pitch salience by assessing the click train tuning
% DEPENDENCIES: 
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, April 2023

function sensitivity = estimate_pitch_salience_sensitivity(CT_tuning,binary_flag,binary_threshold)

% INPUTS:
% CT_tuning: 5 x 17 matrix, where each row is the tuning curve for one of 
%       the CT stimuli IN ORDER so that CT_tuning(1,:) = tuning for CT0,
%       CT_tuning(2,:) = tuning for CT5, etc up to CT40
% binary_flag = 1 or 0, 1 will return the sensitivity as either yes (1)
%       or no (0), 0 will return a continuous metric (optional)
% binary_threshold = floating point to use as the threshold for the slope 
%       if we are returning a binary sensitivity estimate (optional)

if ~exist('binary_flag','var') % if this wasn't passed, assume a continuous returned value
    binary_threshold = 0;
end

if ~exist('binary_threshold','var') % this wasn't passed, assume a middle of the road threshold
    binary_threshold = 0;
end

CT_tuning = zscore(CT_tuning,0,'all'); % zscore the CT_tuning so that the metric is comparable across many neurons

diffs = zeros(size(CT_tuning,1)-1,1); % allocate space to save the difference between the curves
CT0 = CT_tuning(1,:); % extract just the tuning curve for 0 jitter

[~,I] = max(CT0); % find the best frequency (pitch evoking the maximum spike rate)
window = I-4:I+4; % look in the window surrounding the best frequency 

% eliminate pitches outside our stim protocol
window(window<1) = []; 
window(window>17) = [];

% for each CT stimulus, compare the area under the curve to that of CT0
for ct = 1:size(CT_tuning,1)

    diffs(ct) = trapz(CT0(window)) - trapz(CT_tuning(ct,window));
    
end

% fit a line through these points
p = polyfit(1:4,diffs(2:5),1);

if binary_flag % if we want binary, apply the threshold and return
    if p(1) > threshold
        sensitivity = 1;
    else
        sensitivity = 0;
    end
else % else, directly return the slope of the line as the metric
    sensitivity = p(1);
end
end