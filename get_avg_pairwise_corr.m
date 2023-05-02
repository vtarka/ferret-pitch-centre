%% Function to evaluate the average tuning of all the pairs of tuning curves in a response profile
% DEPENDENCIES: 
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, April 2023


function c = get_avg_pairwise_corr(tuning)

% INPUTS:
% tuning - nStims x 17 matrix where each row is the z-scored averaged tuning curve for that stimulus

nStims = size(tuning,1);
corrs = zeros(nchoosek(nStims,2),1);

corr_counter = 1;
for i = 1:nStims
    for k = 1:nStims

        if i<k
            corrs(corr_counter) = corr(tuning(i,:)',tuning(k,:)','rows','complete');
            corr_counter = corr_counter + 1;
        end
    end
end

c = mean(corrs,'omitnan');
end