%% Replicating the key figures from Bendor and Wang, 2005 (Nature) to see if we can get similar results
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, November 2022


%% FINDING THE PITCH NEURONS
%First we find the CF of all neurons
% Then for each neuron, we look at the responses to allHarm (equivalent of B&W's MF) at F0's around the CF
% The allHarm doesn't have to elicit the same CF, only has to elicit a significant response (t-test, p<0.05)

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

neighbor_width = 0; % set how many neighboring frequencies to the BF to incorporate

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

totalPN_count = 0;
PN_CFs = [];
PN_units = cell(length(Animals),4); % allocate space to save the units we find as pitch neurons

% for each recording
for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

    window = [0 0.15]; % response window
    windowSilence = [0.35 0.5]; % inter-trial window

    tone_CF = BFs(:,13); % get the tone-evoked BFs for each unit
    PN_unit_list = [];
    localMFCF = [];

    % for each unit
    for uu = 1:length(units)
        unit = units(uu);
        unitSpikes = Y(Y(:,3)==unit,:); % get the spikes of just this unit

        unitCF = tone_CF(uu); % this unit's tone-evoked BF

        if unitCF==0 % if the unit was F0-sensitive to tones, skip it
            continue
        end

        neighboring_Fs = unitCF-neighbor_width:unitCF+neighbor_width;
        neighboring_Fs(neighboring_Fs<1 | neighboring_Fs>17) = []; % eliminate Fs that weren't presented
        neighboring_Fs = Flist(neighboring_Fs);

        % eliminate neighboring F0s with harmonics less than one octave from unit CF
        for ff = 1:length(neighboring_Fs)
            if neighboring_Fs(ff)*2 < unitCF
                neighboring_Fs(ff) = [];
            end
        end


        % allocate space for spikes in response and silence windows
        nSpikesEvoked = zeros(length(repeats),length(neighboring_Fs));
        nSpikesSilence = zeros(length(repeats),length(neighboring_Fs));
    
        meanSpikes = zeros(length(neighboring_Fs),1);

        % for each frequency
        for ff = 1:length(neighboring_Fs)

            stimNum = find(strcmp(type,'allHarm') & (F0==neighboring_Fs(ff))); % unique name for combination of stim type and F0

            % for each presentation of this stimulus    
            for rr = 1:length(repeats)
    
                spikeIDXs_evoked = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                nSpikesEvoked(rr,ff) = sum(spikeIDXs_evoked);
    
                spikeIDXs_silence = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>windowSilence(1) & unitSpikes(:,2)<windowSilence(2);
                nSpikesSilence(rr,ff) = sum(spikeIDXs_silence);
    
            end % ends repeat loop

            meanSpikes(ff) = mean(nSpikesEvoked(:,ff));

        end % ends F0 loop


        % t-test to see if the neuron responded significantly to allHarm
        h = ttest(nSpikesEvoked(:), nSpikesSilence(:),'Alpha',0.05,'Tail','right');

        % if it passed the t-test, save this unit's info
        if h==1
            totalPN_count = totalPN_count + 1;
            PN_CFs = [PN_CFs, unitCF];
            PN_unit_list = [PN_unit_list, unit];

            [~,i] = max(meanSpikes);
            best_F0 = neighboring_Fs(i);
            F0_IDX = find(Flist==best_F0);
            localMFCF = [localMFCF, F0_IDX];
        end

    end % ends unit loop

    % save all the units we found for this penetration
    PN_units{ap,1} = Animals{ap};
    PN_units{ap,2} = Pens{ap};
    PN_units{ap,3} = PN_unit_list;
    PN_units{ap,4} = localMFCF;

end % ends recording loop


%% Plot PN CFs and non-PN CFs in histogram (BW Fig 2b)

PN_CFs = [];
nPN_CFs = [];
pensOI = 1:20; % penetrations to include in this figure (where each value is a row in PN_units)

% for each penetration
for ap = 1:length(pensOI)
    
    pen = pensOI(ap);
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' PN_units{pen,1} '/tmp/Spikes_' PN_units{pen,1} '_' PN_units{pen,2} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    PNUnits = PN_units{pen,3}; % array of units we found to be pitch neurons

    [~,PNUnit_IDXs] = ismember(PNUnits,allUnits); % find the indices of these pitch neurons within the list of all units

    PN_tone_CFs = BFs(PNUnit_IDXs,13); % get the PN's best frequencies
    PN_CFs = [PN_CFs; PN_tone_CFs];

    nPN_tone_CFs = BFs(:,13); % get the non-PN's best frequencies
    nPN_tone_CFs(PNUnit_IDXs) = [];
    nPN_tone_CFs(nPN_tone_CFs==0) = []; % eliminate any non-PN's that were not F0-sensitive to tones
    nPN_CFs = [nPN_CFs; nPN_tone_CFs];

end % ends recording loop

figure;
[PN_N,PN_edges] = histcounts(PN_CFs,1:18);
PN_CFs_norm = PN_N/length(PN_CFs);

[nPN_N,~] = histcounts(nPN_CFs,1:18);
nPN_CFs_norm = nPN_N/length(nPN_CFs);
colors = colormap(jet(2));
histogram('Categories',string(Flist),'BinCounts',PN_CFs_norm,'FaceColor',colors(1,:),'FaceAlpha',0.5)
hold on
histogram('Categories',string(Flist),'BinCounts',nPN_CFs_norm,'FaceColor',colors(2,:),'FaceAlpha',0.5)
ylim([0 0.30]); yticks(0:0.05:0.3); yticklabels(0:0.05*100:0.3*100)
xlabel('Pure tone CF (Hz)')
ylabel('Percent of samples')
legend({'Pitch Neurons','Non-Pitch Neurons'},'Location','northwest')
set(gca,'FontSize',20)



%% Determine location liklihood of PNs (BW Fig 2a)

penUnits = [27,146,156,211,36,46,9,43,...
    34,77,146,58,60,61,37,29,...
    43,40,26,26]; % the numer of total units in each penetration (in each row of PN_units)

% Make a list of how many PNs per penetration there were
PNUnits = [];

for ap = 1:length(PN_units)
    PNUnits = [PNUnits; length(PN_units{ap,3})];
end

figure;
PNUnitFrequency = PNUnits ./ sum(PNUnits); % percent of units in each penetration determined to be PN_units
plot(PNUnitFrequency,'Color','k','LineWidth',3)
xticks(1:20);xticklabels({'N1','N2','N3','N4','N5','N6','N7','N8','R4','R5','R8','R13','De2','De3','De5','De8','Do0','Do1','Do2','Do4'})
ylim([0 0.35]); yticks(0:0.05:0.35); yticklabels(0:0.05*100:0.35*100)
ylabel('Percent of neurons found to be PNs')
xlabel('Penetration')
set(gca,'FontSize',18)
xlim([1 20])



%% Plot pure tone CFs vs MF CFs (BW Fig 3b)

figure;
CFs = [];

pensOI = [3,4,8,11]; % penetrations to include in this figure (where each value is a row in PN_units)

for ap = 1:length(pensOI)
    
    pen = pensOI(ap);
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' PN_units{pen,1} '/tmp/Spikes_' PN_units{pen,1} '_' PN_units{pen,2} '_Good_Pitch.mat']);

    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    PNUnits = PN_units{pen,3};
    window = [0 0.15];

    PN_MF_CFs = zeros(length(PNUnits),1); % allocate space for MF-evoked BFs

    % for each pitch neuron
    for uu = 1:length(PNUnits)

        unit = PNUnits(uu);
        unitSpikes = Y(Y(:,3)==unit,:); % get the spikes of just this unit

        nSpikesEvoked = zeros(length(repeats),length(Flist));
            
        meanSpikes = zeros(length(Flist),1);

        % for each F0
        for ff = 1:length(Flist)

            stimNum = find(strcmp(type,'allHarm') & (F0==Flist(ff)));  % unique name for combination of stim type and F0

            if isempty(stimNum) % if this sound/F0 combo wasn't presented
                continue
            end
    
            % for each presentation of this stimulus
            for rr = 1:length(repeats)
    
                spikeIDXs_evoked = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                nSpikesEvoked(rr,ff) = sum(spikeIDXs_evoked);
    
            end % ends repeat loop

            meanSpikes(ff) = mean(nSpikesEvoked(:,ff));

        end % ends F0 loop

        [~,i] = max(meanSpikes); % find which F0 evoked the max spike rate
        PN_MF_CFs(uu) = i;
    end

    [~,PNUnit_IDXs] = ismember(PNUnits,allUnits);
    PN_tone_CFs = BFs(PNUnit_IDXs,13);

    CFs = [CFs; PN_tone_CFs, PN_MF_CFs];

end % ends recording loop


figure; hold on
for row = 1:length(CFs)
    x = CFs(row,1);
    y = CFs(row,2);

    repeats = length(find(CFs(:,1)==x & CFs(:,2)==y));
    scatter(x,y, repeats*60,'k','filled')
end

plot([1 17],[1 17],'k--')
xlim([0 18]); xticks (1:17); xticklabels(Flist)
ylim([0 18]); yticks(1:17); yticklabels(Flist)
xlabel('Tone-evoked CF (Hz)')
ylabel('allHarm-evoked CF (Hz)')
text(2,17,'r = 0.37, p = 0.08','FontSize',24)
text(2,16,'n = 24','FontSize',24)
set(gca,'FontSize',20)


%% Plot the normalized discharge rate for every PN over increasing jitter (BW Fig 4a)

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

CT_stims = {'CT0','CT5','CT10','CT20','CT40'}; % stims with jitter
discharge_rates = zeros(32,length(CT_stims),13);
norm_discharge_rates = zeros(32,length(CT_stims));
unit_counter = 1;
neighbor_width = 0; % set how many neighboring frequencies to the BF to incorporate

for ap = 1:length(PN_units)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' PN_units{ap,1} '/tmp/Spikes_' PN_units{ap,1} '_' PN_units{ap,2} '_Good_Pitch.mat']);

    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));

    PNUnits = PN_units{ap,3};
    [~,PNUnit_IDXs] = ismember(PNUnits,allUnits); % find where the pitch neurons lie in the list of all neurons
    tone_CFs = BFs(PNUnit_IDXs,13); % get the pitch neuron's tone-evoked BFs

    window = [0 0.15]; % response window in seconds

    % for each pitch neuron
    for uu = 1:length(PNUnits)

        unit = PNUnits(uu);
        unitSpikes = Y(Y(:,3)==unit,:); % get the spikes of just this unit

        unitCF = tone_CFs(uu); % get this unit's tone-evoked BF

        if unitCF==0 % if this neuron wasn't F0-sensitive in response to tones, skip it
            continue
        end

        neighboring_Fs = unitCF-neighbor_width:unitCF+neighbor_width;
        neighboring_Fs(neighboring_Fs<1 | neighboring_Fs>17) = []; % eliminate Fs that weren't presented
        neighboring_Fs = Flist(neighboring_Fs);


        maxDischarge = 0;
        for ss = 1:length(CT_stims) % for every click-train stimulus

            nSpikes = zeros(length(repeats),length(neighboring_Fs));

            % for each frequency
            for ff = 1:length(neighboring_Fs)
    
                stimNum = find(strcmp(type,CT_stims{ss}) & (F0==neighboring_Fs(ff))); 
              
                for rr = 1:length(repeats)
        
                    spikeIDXs_evoked = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                    nSpikes(rr,ff) = sum(spikeIDXs_evoked);

                    if sum(spikeIDXs_evoked)>maxDischarge % if this discharge rate is the largest we've found so far, save it
                        maxDischarge = sum(spikeIDXs_evoked);
                    end
        
                end % ends repeat loop

            end % ends F0 loop

            discharge_rates(unit_counter,ss,:) = nSpikes;

        end % ends CT stim loop

        discharge_rates(unit_counter,:,:) = discharge_rates(unit_counter,:,:);
        norm_discharge_rates(unit_counter,:) = mean(discharge_rates(unit_counter,:,:),3)./ max(mean(discharge_rates(unit_counter,:,:),3));

        unit_counter = unit_counter + 1; 
    end % ends loop through units
end % ends loop through recordings

figure; hold on
for ct = 1:length(CT_stims)
    for pn = 1:length(discharge_rates)
        scatter(ct,norm_discharge_rates(pn,ct),'k','filled')
    end
end

discharge_rates = mean(discharge_rates,3);
dr_diffs = abs(discharge_rates(:,1)-discharge_rates(:,5));
[~,I] = sort(dr_diffs,'descend');
% I = I(1:13); % uncomment to 'cheat' the figure as B&W did

figure; errorbar(mean(norm_discharge_rates,'omitnan'),std(norm_discharge_rates,'omitnan')/sqrt(length(norm_discharge_rates)),'LineWidth',2,'Color','k')
% ylim([0.05 0.15]); ylabel('Normalized discharge rate')
xlim([0 6]); xticks(1:5); xticklabels([0 5 10 20 40]); xlabel('Maximum jitter (%)')
set(gca,'FontSize',20)

figure; hold on

for ct = 1:length(CT_stims)
    for pn = 1:length(discharge_rates)
        scatter(ct,discharge_rates(pn,ct),'k','filled')
    end
end


%% Discharge rates per lowerst harmonic present (in Hz) (BW Fig 4D)


discharge_rates = zeros(66,2);
allHarmLowest = repmat(2,1,15);
lowLowest = repmat(2,1,15);
highLowest = [5,5,6,7,7,8,8,8,9,9,9,9,8,7,6];

lowestHarm = [allHarmLowest; lowLowest; highLowest];

unit_counter = 1;
skipped_counter = 0;
for ap = 1:length(PN_units)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' PN_units{ap,1} '/tmp/Spikes_' PN_units{ap,1} '_' PN_units{ap,2} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));
    PNUnits = PN_units{ap,3};
    [~,PNUnit_IDXs] = ismember(PNUnits,allUnits);
    tone_CFs = BFs(PNUnit_IDXs,13);

    window = [0 0.15];
    CT_stims = {'allHarm','low','high'};
    tone_CF = BFs(:,13); 

    for uu = 1:length(PNUnits)
        unit = PNUnits(uu);
        unitSpikes = Y(Y(:,3)==unit,:); % get the spikes of just this unit

        unitCF = tone_CFs(uu);
        if (unitCF==1 || unitCF == 2)
            skipped_counter = skipped_counter + 1;
            continue
        end

%         neighboring_Fs = unitCF-0:unitCF+0;
%         neighboring_Fs(neighboring_Fs<1 | neighboring_Fs>17) = []; % eliminate Fs that weren't presented
%         neighboring_Fs = Flist(neighboring_Fs);

        unitDischarge = zeros(length(CT_stims),length(repeats));
        maxDischarge = 0;

        for ss = 1:length(CT_stims)

            nSpikes = zeros(length(repeats),1);

            stimNum = find(strcmp(type,CT_stims{ss}) & (F0==Flist(unitCF))); 

            for rr = 1:length(repeats)
    
                spikeIDXs_evoked = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                nSpikes(rr) = sum(spikeIDXs_evoked);
    
                if sum(spikeIDXs_evoked)>maxDischarge
                    maxDischarge = sum(spikeIDXs_evoked);
                end

            end % ends repeat loop

            unitDischarge(ss,:) = nSpikes;
            discharge_rates(unit_counter,1) = lowestHarm(ss,unitCF-2);
            unit_counter = unit_counter + 1;

        end % ends CT stim loop#

        unitDischarge = unitDischarge; % ./ maxDischarge;
        discharge_rates(unit_counter-ss:unit_counter-1,2) = mean(unitDischarge,2);


    end % ends loop through units
end % ends loop through recordings


fs = unique(discharge_rates(:,1));
DRs = zeros(length(fs),1);
DRs_err = zeros(length(fs),1);
for ff = 1:length(fs)
    f = fs(ff);

    fIDX = find(discharge_rates(:,1)==f);
    DRs(ff) = mean(discharge_rates(fIDX,2),'omitnan');
    DRs_err(ff) = std(discharge_rates(fIDX,2),'omitnan')/sqrt(length(discharge_rates(fIDX,2)));

end

figure;
errorbar(fs,DRs,DRs_err,'Linewidth',2)
% ylim([0 0.45]); 
% xlim([0 27]); xticks(1:2:length(fs)); xticklabels(fs(1:2:end))
% ylabel('Normalized discharge rate')
% xlabel('Frequency of lowest harmonic present (Hz)')
% set(gca,'FontSize',20)



%% Discharge per lowest harmonic present (integer) (BW Fig 4C)

discharge_rates = zeros(66,2);
allHarmLowest = repmat(2,1,17);
lowLowest = repmat(2,1,15);
highLowest = [4,4,5,5,6,7,7,8,8,8,9,9,9,9,8,7,6];

figure; hold on

for ap = 1:length(PN_units)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' PN_units{ap,1} '/tmp/Spikes_' PN_units{ap,1} '_' PN_units{ap,2} '_Good_Pitch.mat']);

    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    allUnits = unique(Y(:,3));
    PNUnits = PN_units{ap,3};
    [~,PNUnit_IDXs] = ismember(PNUnits,allUnits);
    tone_CFs = BFs(PNUnit_IDXs,13);

    window = [0 0.15];
    CT_stims = {'allHarm','low','high'};

    highDischarge = zeros(length(PNUnits),2);
    allHarmDischarge = zeros(length(PNUnits),2);
    lowHarmDischarge = zeros(length(PNUnits),2);

    for uu = 1:length(PNUnits)
        unit = PNUnits(uu);
        unitSpikes = Y(Y(:,3)==unit,:); % get the spikes of just this unit

        unitCF = tone_CFs(uu);

        maxDischarge = 0;

        % HIGH HARM
        nSpikes = zeros(length(repeats),1);

        stimNum = find(strcmp(type,'high') & (F0==Flist(unitCF))); 

        for rr = 1:length(repeats)

            spikeIDXs_evoked = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
            nSpikes(rr) = sum(spikeIDXs_evoked);

            if sum(spikeIDXs_evoked)>maxDischarge
                maxDischarge = sum(spikeIDXs_evoked);
            end

        end % ends repeat loop

        highDischarge(uu,1) = mean(nSpikes);
        highDischarge(uu,2) = highLowest(unitCF);

        % ALL HARM
        nSpikes = zeros(length(repeats),1);

        stimNum = find(strcmp(type,'allHarm') & (F0==Flist(unitCF))); 

        for rr = 1:length(repeats)

            spikeIDXs_evoked = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
            nSpikes(rr) = sum(spikeIDXs_evoked);

            if sum(spikeIDXs_evoked)>maxDischarge
                maxDischarge = sum(spikeIDXs_evoked);
            end

        end % ends repeat loop

        allHarmDischarge(uu,1) = mean(nSpikes);
        allHarmDischarge(uu,2) = allHarmLowest(unitCF);


        if (unitCF==1 || unitCF==2)
            skipped_counter = skipped_counter + 1;
            allHarmDischarge(uu,1) = allHarmDischarge(uu,1) ./ maxDischarge;
            highDischarge(uu,1) = highDischarge(uu,1) ./ maxDischarge;
            lowDischarge(uu,1) = nan;
            continue
        end

        % LOW
        nSpikes = zeros(length(repeats),1);

        stimNum = find(strcmp(type,'low') & (F0==Flist(unitCF))); 

        for rr = 1:length(repeats)

            spikeIDXs_evoked = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
            nSpikes(rr) = sum(spikeIDXs_evoked);

            if sum(spikeIDXs_evoked)>maxDischarge
                maxDischarge = sum(spikeIDXs_evoked);
            end

        end % ends repeat loop

        lowDischarge(uu,1) = mean(nSpikes);
        lowDischarge(uu,2) = lowLowest(unitCF-2);

        allHarmDischarge(uu,1) = allHarmDischarge(uu,1) ./ maxDischarge;
        highDischarge(uu,1) = highDischarge(uu,1) ./ maxDischarge;
        lowDischarge(uu,1) = lowDischarge(uu,1) ./ maxDischarge;

    end % ends loop through units

    for uu = 1:length(PNUnits)
        scatter(allHarmDischarge(uu,2),allHarmDischarge(uu,1),'k','filled')
        scatter(lowDischarge(uu,2),lowDischarge(uu,1),'k','filled')
        scatter(highDischarge(uu,2),highDischarge(uu,1),'k','filled')
    end
end % ends loop through recordings


% fs = unique(discharge_rates(:,1));
% DRs = zeros(length(fs),1);
% DRs_err = zeros(length(fs),1);
% for ff = 1:length(fs)
%     f = fs(ff);
% 
%     fIDX = find(discharge_rates(:,1)==f);
%     DRs(ff) = mean(discharge_rates(fIDX,2),'omitnan');
%     DRs_err(ff) = std(discharge_rates(fIDX,2),'omitnan')/sqrt(length(discharge_rates(fIDX,2)));
% 
% end
% 
% figure;
% errorbar(fs,DRs,DRs_err,'Linewidth',2)
ylim([0 0.6]); 
xlim([1 10]); xticks(2:2:10); 
ylabel('Normalized discharge rate')
xlabel('Lowest harmonic present')
set(gca,'FontSize',20)

