%% Function to assess how sensitive a neuron is to phase shifts (alternating and random)
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, April 2023


function sensitivity = estimate_phase_shift_sensitivity(phase_tuning)

% INPUTS:
% phase_tuning: 3 x 17 matrix where each row is the tuning curve for
%        1) high harmonics, 2) alt phase, and 3) rand phase

% Find the lag at which the correlation between no shift and alt shift is highest
[r_alt,lags_alt] = xcorr(phase_tuning(1,:),phase_tuning(2,:));
[~,max_r_idx] = max(r_alt);
best_lag_alt = lags_alt(max_r_idx);

[r_rand,lags_rand] = xcorr(phase_tuning(1,:),phase_tuning(3,:));
[~,max_r_idx] = max(r_rand);
best_lag_rand = lags_rand(max_r_idx);

% Take the z-score of all the tuning
zscored_tuning = zscore(phase_tuning,0,'all');
rand_zscored_tuning = zscored_tuning(3,:);

% if the alt phase tuning is shifted towards lower frequencies, AND
% the rand phase tuning is either flat OR also shifted, it is sensitive
if best_lag_alt > 1 %% && (max(rand_zscored_tuning) < 0.5 || best_lag_rand > 0)
    sensitivity = 1;
else
    sensitivity = 0;
end
end