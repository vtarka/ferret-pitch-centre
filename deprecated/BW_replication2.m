%% First we find the CF of all neurons
% Then for each neuron, we look at the responses to 'low' and 'high' harmonic
% conditions at F0's around the CF
% The MF stimuli don't have to elicit the same CF, only has to elicit a
% significant response (2 SD above spontaneous mean)

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

totalPN_count = 0;
PN_CFs = [];
PN_units = cell(length(Animals),4);

for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    window = [0 0.15];
    windowSilence = [0.35 0.5];

    tone_CF = BFs(:,13);
    PN_unit_list = [];
    localMFCF = [];

    for uu = 1:length(units)
        unit = units(uu);
        unitSpikes = Y(Y(:,3)==unit,:); % get the spikes of just this unit

        unitCF = tone_CF(uu);

        if unitCF==0
            continue
        end

        neighboring_Fs = unitCF-0:unitCF+0;
        neighboring_Fs(neighboring_Fs<1 | neighboring_Fs>17) = []; % eliminate Fs that weren't presented
        neighboring_Fs = Flist(neighboring_Fs);

        nSpikesEvoked = zeros(length(repeats),length(neighboring_Fs));
        nSpikesSilence = zeros(length(repeats),length(neighboring_Fs));
    
        active = 1; % flag to keep track of whether the unit is responsive to both low and high harms

        meanSpikes = zeros(length(neighboring_Fs),1);
        for ff = 1:length(neighboring_Fs)

            if (neighboring_Fs(ff)==250 || neighboring_Fs(ff)==297)
                continue
            end

            stimNum = find(strcmp(type,'low') & (F0==neighboring_Fs(ff))); 
    
            for rr = 1:length(repeats)
    
                spikeIDXs_evoked = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                nSpikesEvoked(rr,ff) = sum(spikeIDXs_evoked);
    
                spikeIDXs_silence = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>windowSilence(1) & unitSpikes(:,2)<windowSilence(2);
                nSpikesSilence(rr,ff) = sum(spikeIDXs_silence);
    
            end % ends repeat loop

            meanSpikes(ff) = mean(nSpikesEvoked(:,ff));

        end % ends F0 loop

        spontaneous_m = mean(nSpikesSilence(:));
        spontaneous_std = std(nSpikesSilence(:));

        evoked_m = mean(nSpikesEvoked(:));

        if evoked_m < spontaneous_m + 2*spontaneous_std
            active = 0;
            continue
        end

        meanSpikes = zeros(length(neighboring_Fs),1);
        nSpikesEvoked = zeros(length(repeats),length(neighboring_Fs));
        nSpikesSilence = zeros(length(repeats),length(neighboring_Fs));
        for ff = 1:length(neighboring_Fs)

            stimNum = find(strcmp(type,'high') & (F0==neighboring_Fs(ff))); 

            if (neighboring_Fs(ff)==250 || neighboring_Fs(ff)==297)
                continue
            end
    
            for rr = 1:length(repeats)
    
                spikeIDXs_evoked = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                nSpikesEvoked(rr,ff) = sum(spikeIDXs_evoked);
    
                spikeIDXs_silence = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>windowSilence(1) & unitSpikes(:,2)<windowSilence(2);
                nSpikesSilence(rr,ff) = sum(spikeIDXs_silence);
    
            end % ends repeat loop

            meanSpikes(ff) = mean(nSpikesEvoked(:,ff));

        end % ends F0 loop

        spontaneous_m = mean(nSpikesSilence(:));
        spontaneous_std = std(nSpikesSilence(:));

        evoked_m = mean(nSpikesEvoked(:));

        if evoked_m < spontaneous_m + 2*spontaneous_std
            active = 0;
            continue
        end

        totalPN_count = totalPN_count + 1;
        PN_CFs = [PN_CFs, unitCF];
        PN_unit_list = [PN_unit_list, unit];

%             [~,i] = max(meanSpikes);
%             best_F0 = neighboring_Fs(i);
%             F0_IDX = find(Flist==best_F0);
%             localMFCF = [localMFCF, F0_IDX];

    end % ends unit loop

    PN_units{ap,1} = Animals{ap};
    PN_units{ap,2} = Pens{ap};
    PN_units{ap,3} = PN_unit_list;

end % ends recording loop

