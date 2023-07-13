

function response = get_binned_signal(unitSpikes, stimulus_number, repeats, bin_size, window)

stimNum = stimulus_number;
nRepeats = length(repeats);

bins = min(window):bin_size:max(window);

nSpikes = zeros(nRepeats,length(bins)-1);
for rr = 1:nRepeats

    for bb = 1:length(bins)-1

        spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>bins(bb) & unitSpikes(:,2)<bins(bb+1);
        nSpikes(rr,ff) = sum(spikeIDXs);

    end % ends bin loop
end % ends repeat loop

response = mean(nSpikes);

end