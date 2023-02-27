%% Counts the number of peaks in a tuning curve
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, January 2023

% 75 (or input threshold) % of maximally evoked activity is a peak
% CASE 1: entries on either side of a peak either don't exist (at edge) or
% are < 75% of the maximum
% CASE 2: entry on left side of a peak is < 75% max (opens parenthesis)
% CASE 3: entry on right side of peak is < 75% max (closes parenthesis)

function [nPeaks, peak_idxs] = count_peaks(tuning_curve, threshold)

max_val = max(tuning_curve); % maximum value in the tuning curve

adjusted_thresh = max_val * threshold; % threshold for another value to be considered a peak

nPeaks = 0;
open = 0; 
open_idx = 0; 
peak_idxs = [];

for i = 1:length(tuning_curve)
    
    % if this response is above threshold
    if tuning_curve(i) > adjusted_thresh 

        % if this is the first value is above threshold
        if i == 1 
            if tuning_curve(i+1) < adjusted_thresh % if the next value is below threshold, the first value is a peak by itself
                peak_idxs = [peak_idxs; i];
                nPeaks = nPeaks + 1;
            else % if the next value is above threshold, the first value is just the start of the peak
                open = 1; % indicate we're currently evaluating an ongoing peak
                open_idx = i;
            end

        % if the last value is above threshold
        elseif i == length(tuning_curve) 
            if tuning_curve(i-1) > adjusted_thresh % if the preceding value is above peak, we are closing an ongoing peak
                nPeaks = nPeaks + 1;
                [~,peak_within] = max(tuning_curve(open_idx:i)); 
                peak_idxs = [peak_idxs; open_idx + peak_within - 1]; % save the maximum value within a multi-index peak as the peak
            else
                peak_idxs = [peak_idxs; i];
                nPeaks = nPeaks + 1;
            end

        % if the values on either side are below threshold, this is a single index peak    
        elseif tuning_curve(i-1) < adjusted_thresh && tuning_curve(i+1) < adjusted_thresh 
            peak_idxs = [peak_idxs; i];
            nPeaks = nPeaks + 1;

        % if only the preceding value is below threshold, this is the start of a multi-index peak
        elseif tuning_curve(i-1) < adjusted_thresh 
            open = 1;
            open_idx = i;

        % if the next value is below threshold, we've reached the end of a multi-index peak        
        elseif tuning_curve(i+1) < adjusted_thresh && open == 1 
            nPeaks = nPeaks + 1;
            [~,peak_within] = max(tuning_curve(open_idx:i));
            peak_idxs = [peak_idxs; open_idx + peak_within - 1]; % save the maximum value within a multi-index peak as the peak
            open = 0;
        end

    end
end