%% Counts the number of peaks in a tuning curve
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, January 2023

% 75% of maximally evoked activity is a peak
% CASE 1: entries on either side of a peak either don't exist (at edge) or
% are < 75% of the maximum
% CASE 2: entry on left side of a peak is < 75% max (opens parenthesis)
% CASE 3: entry on right side of peak is < 75% max (closes parenthesis)

function [nPeaks, peak_idxs] = count_peaks(tuning_curve, threshold)

max_val = max(tuning_curve);

adjusted_thresh = max_val * threshold;

nPeaks = 0;
open = 0;
open_idx = 0;
peak_idxs = [];

%% TODO: implement an indexes method so that the middle of a multi-value peak is returned

for i = 1:length(tuning_curve)
    
    if tuning_curve(i) > adjusted_thresh

        if i == 1
            if tuning_curve(i+1) < adjusted_thresh
                peak_idxs = [peak_idxs; i];
                nPeaks = nPeaks + 1;
            else
                open = 1;
                open_idx = i;
            end
        elseif i == length(tuning_curve)
            peak_idxs = [peak_idxs; i];
            nPeaks = nPeaks + 1;
        elseif tuning_curve(i-1) < adjusted_thresh && tuning_curve(i+1) < adjusted_thresh
            peak_idxs = [peak_idxs; i];
            nPeaks = nPeaks + 1;
        elseif tuning_curve(i-1) < adjusted_thresh
            open = 1;
            open_idx = i;
        elseif tuning_curve(i+1) < adjusted_thresh && open == 1
            nPeaks = nPeaks + 1;
            [~,peak_within] = max(tuning_curve(open_idx:i));
            peak_idxs = [peak_idxs; open_idx + peak_within - 1];
            open = 0;
        end

    end
end
end