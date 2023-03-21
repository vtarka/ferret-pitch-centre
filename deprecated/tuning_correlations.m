

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

for ap = 1:length(PN_units)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' PN_units{ap,1} '/tmp/Spikes_' PN_units{ap,1} '_' PN_units{ap,2} '_Good_Pitch.mat']);

    % function a = plot_tuning_by_cond(Y,type,F0,unit,stims,BFs,animal,pen)

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));
    PNUnits = PN_units{ap,3};
    [~,PNUnit_IDXs] = ismember(PNUnits,allUnits);
    BFs = BFs(PNUnit_IDXs,[13,8,1]);
    window = [0 0.15];

    stims_to_plot = {'CT0','low','high'};

    for uu = 1:length(PNUnits)
        unit = PNUnits(uu);
        unitSpikes = Y(Y(:,3)==unit,:); % get the spikes of just this unit

        tuning = zeros(length(stims_to_plot),length(Flist));
        for ss = 1:length(stims_to_plot)
            nSpikes = zeros(length(repeats),length(Flist)); %creates an array of length repList by length Flist
            
            for ff = 1:length(Flist) % do the below for all the frequencies
        
                stimNum = find(strcmp(type,stims{ss}) & (F0==Flist(ff))); 
    
                if isempty(stimNum) % if this stim type and fo combo wasn't presented
                    nSpikes(:,ff) = 0;
                    continue
                end
                
                % finds the stim label that corresponds to this stim type at
                % this particular F0
    
                % this stimNum will have been presented multiple times, 
                % so go through each presentation
        
                for rr = 1:length(repeats)
                
                    spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                    nSpikes(rr,ff) = sum(spikeIDXs);
    
                end   
            end
    
            nSpikes = mean(nSpikes ./ diff(window)); % spikes per second
            tuning(ss,:) = nSpikes;

        end % end stim loop

        plot_tuning_by_cond(Y,type,F0,unit,stims_to_plot,BFs,PN_units{ap,1},PN_units{ap,2});


        figure('Position',[200 500 1800 1200])
        cc = corrcoef(tuning');
        cc = cc - diag(diag(cc));
        imagesc(cc); colorbar
        xticks(1:5); xticklabels(stims_to_plot)
        yticks(1:5); yticklabels(stims_to_plot)

        sgtitle(sprintf('%s, %s unit # %d',PN_units{ap,1},PN_units{ap,2},unit))

        pause 

        close all

    end % ends unit loop
end % ends recording loop
