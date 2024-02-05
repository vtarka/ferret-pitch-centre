%% Function to assess how sensitive a neuron is to phase shifts (alternating and random)
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, April 2023


function sensitivity = estimate_phase_shift_sensitivity_allHarm(phase_tuning)

% INPUTS:
% phase_tuning: 3 x 17 matrix where each row is the tuning curve for
%        1) high harmonics, 2) alt phase, and 3) rand phase

[~,high_BF] = max(phase_tuning(1,:));
[~,allHarm_BF] = max(phase_tuning(2,:));

% Find the lag at which the correlation between no shift and alt shift is highest
[r_allHarm,lags_allHarm] = xcorr(phase_tuning(1,:),phase_tuning(2,:));
[~,max_r_idx] = max(r_allHarm);
best_lag_allHarm = lags_allHarm(max_r_idx);

% if the alt phase tuning is shifted towards lower frequencies, AND
% the rand phase tuning is either flat OR also shifted, it is sensitive
if best_lag_allHarm > 1 && high_BF~=allHarm_BF %% && (max(rand_zscored_tuning) < 0.5 || best_lag_rand > 0)
    sensitivity = 1;
else
    sensitivity = 0;
end
end