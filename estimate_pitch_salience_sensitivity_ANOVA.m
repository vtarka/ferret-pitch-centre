%% Function to measure how sensitive a neuron is to pitch salience by assessing the click train tuning
% DEPENDENCIES: polyfit
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, April 2023

function sensitivity = estimate_pitch_salience_sensitivity_ANOVA(CT_tuning,alpha)

% INPUTS:
% CT_tuning: 5 x 1 cell array, where CT_tuning{1} holds an nRepeats x
% nFreqs matrix of spikerates in response to CT0

% CT_tuning(3) = [];
CT0 = CT_tuning{1}; % extract just the tuning curve for 0 jitter
CT0_tuning = mean(CT0);

[~,I] = max(CT0_tuning); % find the best frequency (pitch evoking the maximum spike rate)
window = I-2:I+2; % look in the window surrounding the best frequency
window(window<1) = [];
window(window>17) = [];
nFreqs = length(window);

data = zeros(13*nFreqs*5,1);
jitter = zeros(13*nFreqs*5,1);
freq = zeros(13*nFreqs*5,1);

lp = 1;

for ss = 1:length(CT_tuning) % for each jitter level (0 to 40)

    responses = CT_tuning{ss};

    for ff = 1:length(window)

        f_responses = responses(:,window(ff));

        nResponses = length(f_responses);

        data(lp:lp+nResponses-1) = f_responses;
        freq(lp:lp+nResponses-1) = ff;

        if ss<4
            jitter(lp:lp+nResponses-1) = 1;
        else
            jitter(lp:lp+nResponses-1) = 2;
        end

        lp = lp+nResponses;
    end
end

p = anovan(data,{freq,jitter},'model','interaction','display','off');

if p(2) < alpha || p(3) < alpha
    sensitivity = 1;
    return
end

% removing CT10
% CT_tuning(3) = [];

data = zeros(13*nFreqs*5,1);
jitter = zeros(13*nFreqs*5,1);
freq = zeros(13*nFreqs*5,1);

lp = 1;

for ss = 1:length(CT_tuning) % for each jitter level (0 to 40)

    responses = CT_tuning{ss};

    for ff = 1:length(window)

        f_responses = responses(:,window(ff));

        nResponses = length(f_responses);

        data(lp:lp+nResponses-1) = f_responses;
        freq(lp:lp+nResponses-1) = ff;

        if ss<3
            jitter(lp:lp+nResponses-1) = 1;
        else
            jitter(lp:lp+nResponses-1) = 2;
        end

        lp = lp+nResponses;
    end
end

p = anovan(data,{freq,jitter},'model','interaction','display','off');

if p(2) < alpha || p(3) < alpha
    sensitivity = 1;
else
    sensitivity = 0;
end

end