%% Makes all the plots on the FENS 2019 poster

addpath('E:\Veronica\Daniel_ephys_analysis\pitch_ephys')
addpath('E:\Veronica\Daniel_ephys_analysis\pitch_ephys\plotFunc')
addpath(genpath('E:\Veronica\Daniel_ephys_analysis\pitch_ephys\2pi_code_tonotopy\2pi_code_tonotopy'))
addpath(genpath('E:\Veronica\Dan_docs\MATLAB\Project\npy-matlab-master\npy-matlab-master'))

cd('E:\Veronica\Daniel_ephys_analysis\pitch_ephys\DansMATLABData');
load('allDat_allanimalsandpens.mat'); %load all the data with Noah
load('All1wayanova.mat'); %load all the data with Noah

%% Edit spike mat to add stim info. Only need to do this once!

%stimuli are 200ms long for all ferrets apart from Noah - these are 300ms
cd('E:\Veronica\Daniel_ephys_analysis\pitch_ephys') %go to the directory for the code
clear; %get rid of previously loaded stuff
close all; %get rid of diagrams

for g = 1:8 %for each penetration
    
    g
    Animal = 'Noah'; %which animal is it
    Pen = ['P0' num2str(g)]; %which pen is it 
    Qualia = 'Good';
    % Qualia = 'MUA'; %choose which type of neuron activity it is when adding data

    %addpath('D:\Work\Code\gitDev\npy-matlab-master')


    % dataPath = 'D:\Work\Data\ephys\P05-pitch2018\P05-pitch2018_sorted';
    %dataPath = 'D:\Work\Data\ephys\Derry\P02\P02-pitch\CRA';
    % dataPath = '/Volumes/groups/king-archive/Quentin/pitch ephys curated/Derry/P02/P02-pitch/CRA'; %where are the specific data?
    dataPath = ['E:\Noah\For_analysis\P0' num2str(g) '\P0' num2str(g) '-pitch70dB2020_g0'];
    % load('D:\Work\Data\ephys\P05-pitch2018\P05-pitch2018_sorted\gridInfo.mat');
    load([dataPath '\gridInfo.mat']); %loads the file which contains ??
    synch_ch = get_synch(dataPath);  % Get sweep start & stop
    %synch_ch2 = get_synch(dataPath); % Kerry and Dan just testing stuff to
    %see how synch works.

    min_trig_length_s = 0.5; %set to 0.5 for Noah. Was 0.7 for other ferrets
    min_interTrig_length_s = 0.05;
    fs = 30000; %Hz

    disp('Getting triggers')
    [start_time_ms] = get_triggers(synch_ch,min_trig_length_s,min_interTrig_length_s,fs); % Get triggers

    disp('Getting spike times')
    Y = get_spike_times(dataPath,Qualia,start_time_ms);
    %[Y] = loadSpikeTimes(dataPath,Qualia,start_time_ms,fs); %creates an array with functioned version of: CRA, MUA/good, when time started, sampling frequency)
    % Y -> [spike absolute time - spike relative times - unit - stimulus # -
    % repeat # - sweep #] these are the values that keep coming up later
    % Grid info
    stimNames = grid.stimFiles; %gets all the types of stimulus into a variable
    
    r = regexp(stimNames,'_(\d*)k?Hz','tokens'); %makes a variable looking for each stimulus type, an underscore, a number in Hz, and saves any appropriate variables to a token
    F0={};
    for ii = 1:length(r) %for all of the appropriate stimulus files
        if isempty(r{ii}) %if the file is empty
            F0{ii} = '0'; %frequency for this file is 0
        else
            F0{ii} = r{ii}{1}{1}; %otherwise, frequency of this file is equal to the found value in the file
        end
    end
    F0 = F0'; %makes the conjugate complex transpose of the frequency matrix (makes it 1 by 3?)
    F0 = str2double(F0); %makes frequency a number

    type = regexp(stimNames,'_(\w*)_\d*k?Hz|_(\w*)_','tokens'); %looks for file types that are a stimulus name, underscore, any text, a frequency value in Hz, any text, underscore, and makes it a token
    type = cellfun(@(x)((x{1})),type)'; %looks at the first part of the array which in this case is the stimulus name
    type = cellfun(@(x)([x]),type,'UniformOutput',false); %returns the outputs in a cell array with x being the file values that were specified, and y being the stimulus type
    
    for c = 1:length(type) %removing the extra stimuli for Noah - Dan
        if strcmp(type{c},'SAMtones') == 1      
            type{c} = {};
        elseif strcmp(type{c},'allHarmRand') == 1
            type{c} = {};
        else
            type{c} = type{c};
        end
    end
    
    type = type(~cellfun(@isempty,type));

    disp('Saving file')
%     save(['Spikes_' Animal '_' Pen '_' Qualia '_Pitch.mat'],'Y','type','F0','dataPath','start_time_ms'); %save variables that have been created by the above code
    
end

%% Remove extra frequencies from Noah

for g = 1:8 %for each penetration
    
    g
    Animal = 'Noah'; %which animal is it
    Pen = ['P0' num2str(g)]; %which pen is it 
    Qualia = 'Good';
    cd(['E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData\' Animal])
    load(['E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData\' Animal '\Spikes_' Animal '_' Pen '_' Qualia '_Pitch.mat']);

    F0(118:151) = []; %frequencies for SAMtones and allHarmRand
%     save(['Spikes_' Animal '_' Pen '_' Qualia '_Pitch.mat'],'Y','type','F0','dataPath','start_time_ms'); %save variables that have been created by the above code

end

%% Plot responses to Pitch for 1 stim through units

clear; %resets any previous values
close all %closes all previously open figures whose handles can be seen

Animal = 'Noah';
Pen = 'P02';
Qualia = 'Good';
load(['E:\Veronica\Daniel_ephys_analysis\pitch_ephys\DansMATLABData\' Animal '\Spikes_' Animal '_' Pen '_' Qualia '_Pitch.mat']);

stim = 'allHarm'; %selects stimulus to be plotted
Flist = unique(F0); %creates a list of F0 values without repetitions
UnitList = unique(Y(:,3)); %creates a list of non-repeating units
repList = unique(Y(:,5)); %creates a list of individual trials
window = [0 0.5]; %seconds %why this long? Equivalent for Noah?

a1 = figure; %this just creates a figure object so a specific one can be specified later
a2 = figure;
% disp(length(UnitList))
p=1;
total = 0;
for i = 1:length(UnitList), %do the below between the value of 9 and the last unique unit
% for i = 15:157
    Nspikes = zeros(length(repList),length(Flist)); %creates an array of length repList by length Flist
    SpikeTimes = {}; %blank array 
    c = 1; %will be useful later
    for f = 1:length(Flist), %do the below for all the frequencies
        for r = 1:length(repList), %do the below for all the spikes
            stimNum = find(strcmp(type,stim) & (F0==Flist(f))); %creates an array with the position of the correct stimulus and the currently investigated frequency
            idx = (Y(:,4) == stimNum) & (Y(:,3) == UnitList(i)) & (Y(:,5)==repList(r)) & (Y(:,2) > window(1)) & (Y(:,2) < window(2)); %creates a valid spike variable if correct conditions are satisfied
            Nspikes(r,f) = sum(idx);
            SpikeTimes{c} = Y(idx,2);
            c = c+1;
        end
    end
    Nspikes = Nspikes ./ diff(window); %Spikes per second
    spikes = sum(Nspikes,'all');
    meanspikes = spikes/(length(Flist)*length(repList));
    total = total + meanspikes; %divides the number of spikes by the time window to make them a rate
    figure(a1); %the following is for figure a1
    clf %clear previous figure
%     subplot(3,3,p)
    errorbar(Flist,mean(Nspikes),ste(Nspikes)) %error bar plot for X axis frequency, y axis mean spike number, error bars are standard error for the mean spike number
    title(sprintf('Unit number %d',UnitList(i))) %title is this for each unit number
%     txt = ['n = ' (num2str(Nspikes(r)*(Nspikes(f))))];
%     text(1000,1, txt)
    set(gca,'XScale','log') %make the x axis logarithmic
    xlabel('F0 (Hz)') %x axis label is the F0
    ylabel('Evoked Firing Rate (spike/sec)'); %y axis label is the firing rate
    
    figure(a2); %the following is for figure a2
    clf %clear previous a2
%     subplot(3,3,p)
    hold on %plot this stuff and don't lose it
    for k = 1:length(SpikeTimes) %do the following for all the spike times    
        s = scatter(SpikeTimes{k},ones(length(SpikeTimes{k}),1) * k,'marker','o','markerfacecolor','k','markeredgecolor','k'); %create a scatter plot of spike times versus
        s.SizeData = 10; %make the dots in the raster plot that size
    end
    y = length(repList):length(repList):length(repList)*(length(Flist)); %Y is a 3D vector where X and Y equal ? and Z is this * number of frequencies
    line(repmat(window,length(y),1)',[y' y']','color','r');
    hold off
    title(sprintf('Unit number %d',UnitList(i)))
    xlabel('Time relative to stimulus onset (second)')
    ylabel('F0 (Hz)')
    set(gca,'YTick',y)
    set(gca,'YTickLabel',Flist)
    p=p+1;
    pause();
end
overallmean = total/length(Nspikes);
disp(overallmean)

%% Plot responses to Pitch for 1 units through all stim

clear;
close all
Animal = 'Noah';
Pen = 'P02';
Qualia = 'Good';
load(['E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData\' Animal '\Spikes_' Animal '_' Pen '_' Qualia '_Pitch.mat']);


Unit = 419;
stimList = unique(type); %creates a list of all the unique stimuli
Flist = unique(F0);
UnitList = unique(Y(:,3));
repList = unique(Y(:,5));
window = [0 0.5]; %seconds

b1 = figure;
b2 = figure;



for i = 1:length(stimList),
    stim = stimList{i};
    Nspikes = zeros(length(repList),length(Flist));
    SpikeTimes = {};
    c = 1;
    for f = 1:length(Flist),
        for r = 1:length(repList),
            stimNum = find(strcmp(type,stim) & (F0==Flist(f)));
            if isempty(stimNum),
                fprintf('Stim : %s , F0: %d',stim,Flist(f));
                continue;
            end
            idx = (Y(:,4) == stimNum) & (Y(:,3) == Unit) & (Y(:,5)==repList(r)) & (Y(:,2) > window(1)) & (Y(:,2) < window(2));
            Nspikes(r,f) = sum(idx);
            SpikeTimes{c} = Y(idx,2);
            c = c+1;
        end
    end

    Nspikes = Nspikes ./ diff(window);
    figure(b1);
    clf
    errorbar(Flist,mean(Nspikes),ste(Nspikes))
    title(sprintf('Unit: %d - %s',Unit,stim))
    set(gca,'XScale','log')
    set(gcf,'position',[50 500 477 339])
    xlabel('F0 (Hz)')
    ylabel('Evoked Firing Rate (spike/sec)');
    %saveas(a1,sprintf('./Figures2/Tuning/U%d_%s_Selectivity.png',Unit,stim));
    
    figure(b2);
    clf
    hold on
    for k = 1:length(SpikeTimes)
        s = scatter(SpikeTimes{k},ones(length(SpikeTimes{k}),1) * k,'marker','o','markerfacecolor','k','markeredgecolor','k');
        s.SizeData = 10;
    end
    y = length(repList):length(repList):length(repList)*(length(Flist));
    line(repmat(window,length(y),1)',[y' y']','color','r');
    hold off
    title(sprintf('Unit: %d - %s',Unit,stim))
    set(gcf,'position',[600 500 477 339])
    xlabel('Time relative to stimulus onset (second)')
    ylabel('F0 (Hz)')
    set(gca,'YTick',y)
    set(gca,'YTickLabel',Flist)
    %saveas(b2,sprintf('./Figures2/Raster/U%d_%s_Raster.png',Unit,stim));
    
     pause();
end

%% Save Tuning curves

close all

Animal = 'Derry';
Pen = 'P02';
Qualia = 'Good';
load(['E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData\' Animal '\Spikes_' Animal '_' Pen '_' Qualia '_Pitch.mat']);

%load('Spikes_Derry_P02_Pitch.mat.mat');
%load('E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\Quentin Pitch ephys\PitchEphysData\Spikes_Derry_P02_Pitch.mat');

% stim = 'allHarm';
% stim2 = 'tone';
stimList = unique(type); %makes a list of all the unique stimulus types
Flist = unique(F0); %makes a list of all the unique fundamental frequencies
UnitList = unique(Y(:,3));
repList = unique(Y(:,5));
window = [0 0.15]; %seconds

c1 = figure;
% h2 = figure;

% allResp = nan(length(UnitList),length(stimList)+1);

for i = 1:length(UnitList),
    subY = Y(Y(:,3) == UnitList(i),:);
%     allResp(i,1) = UnitList(i);
    for s = 1:length(stimList),
        stim = stimList{s};
        Nspikes = zeros(length(repList),length(Flist));
        SpikeTimes = {};
        c = 1;
        for f = 1:length(Flist),
            for r = 1:length(repList),
                stimNum = find(strcmp(type,stim) & (F0==Flist(f)));
                if isempty(stimNum), continue; end;
%                 idx = (Y(:,4) == stimNum) & (Y(:,3) == UnitList(i)) & (Y(:,5)==repList(r)) & (Y(:,2) > window(1)) & (Y(:,2) < window(2));
                idx = (subY(:,4) == stimNum) & (subY(:,5)==repList(r)) & (subY(:,2) > window(1)) & (subY(:,2) < window(2));
                Nspikes(r,f) = sum(idx);
                SpikeTimes{c} = subY(idx,2);
                c = c+1;
            end
        end
        Nspikes= Nspikes./ diff(window);
        
        [p,tbl,stats] = anova1(Nspikes,[],'off');
        
        figure(c1);
        clf
        errorbar(Flist,mean(Nspikes),ste(Nspikes))
        set(gca,'XScale','log')
        xlabel('F0 (Hz)')
        ylabel('Evoked Firing Rate (spike/sec)');
        title(sprintf('Unit %d | p=%0.2f | %s',UnitList(i),p,stim))
%         saveas(c1,sprintf('./Figures2/TuningCurves/%s_%s_U%d_%s_TC.png',Animal,Pen,UnitList(i),stim));
        pause();
    end
    
end

%% T test for responsiveness
% ANOVA for all stimuli
% Save all figures + resp matrix
% cd('E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData')

clear;
close all

% Animal = 'Ronnie';
% Pen = 'P13';
% 
% load(['Spikes_' Animal '_' Pen '_Pitch.mat']);
%for g = 1:8
for g = 13 
   
    Animal = 'Ronnie';
    Pen = ['P' num2str(g)];
    Qualia = 'Good';
    load(['E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData\' Animal '\Spikes_' Animal '_' Pen '_' Qualia '_Pitch.mat']);

    cd (['E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData\' Animal])

%load('E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData\tempBin\Spikes_Noah_P01_Good_Pitch.mat');
    % stim = 'allHarm';
    % stim2 = 'tone';
    stimList = unique(type);
    Flist = unique(F0);
    UnitList = unique(Y(:,3));
    repList = unique(Y(:,5));
    window = [0 0.15]; %seconds
    windowSilence = [0.3 0.4]; % seconds %Dan changed from 0.8-0.95 to 0.4-0.5 for Noah to reflect timing differences 
                               %this will need to change for other ferrets
                               %too, to be equivalent (0.3-0.4)
                               
    % allResp = nan(length(UnitList),length(stimList)+1);

    for i = 1:length(UnitList),
        subY = Y(Y(:,3) == UnitList(i),:);
        allNspikes = [];
        allNspikesSilence = [];
        pFullBreak = [];

        for s = 1:length(stimList),
            stim = stimList{s};
            Nspikes = zeros(length(repList),length(Flist));
            NspikesSilence = zeros(length(repList),length(Flist));
            c = 1;
            for f = 1:length(Flist), % 
                for r = 1:length(repList),
                    stimNum = find(strcmp(type,stim) & (F0==Flist(f)));
                    if isempty(stimNum), disp(num2str(stimNum)); continue; end;
                    idx = (subY(:,4) == stimNum) & (subY(:,5)==repList(r)) & (subY(:,2) > window(1)) & (subY(:,2) < window(2));
                    Nspikes(r,f) = sum(idx);
                    idx = (subY(:,4) == stimNum) & (subY(:,5)==repList(r)) & (subY(:,2) > windowSilence(1)) & (subY(:,2) < windowSilence(2));
                    NspikesSilence(r,f) = sum(idx);
                    c = c+1;
                end
            end 
            Nspikes= Nspikes./ diff(window);
            NspikesSilence = NspikesSilence./ diff(windowSilence);

            [~,p] = ttest(Nspikes,NspikesSilence,0.05,'right');
            pFullBreak = [pFullBreak p];

            [~,p] = ttest(Nspikes(:),NspikesSilence(:),0.05,'right');
            allResp(i).PerStim.p(s) = p;
            allResp(i).PerStim.N = length(Nspikes(:));

            allNspikes = [allNspikes; Nspikes];
            allNspikesSilence = [allNspikesSilence; NspikesSilence];
        end
        [~,p] = ttest(allNspikes(:),allNspikesSilence(:),0.05,'right');
        allResp(i).Global.p = p;
        allResp(i).Global.N = length(allNspikes(:));

        [~,p] = ttest(allNspikes,allNspikesSilence,0.05,'right'); 
        allResp(i).PerFreq.p = p;
        allResp(i).PerFreq.N = length(allNspikes);

        pFullBreak(isnan(pFullBreak)) = 1;
        allResp(i).FullBreak.p = pFullBreak;
        allResp(i).FullBreak.N = length(repList);
    end

    save(['Responsiveness_' Animal '_' Pen '_Pitch.mat'],'allResp','stimList','Flist','UnitList');
    disp('Done.')
end

%% Two-way ANOVA for Selectivity

clear;
close all

Animal = {'Ronnie' 'Ronnie' 'Ronnie' 'Ronnie' 'Derry' 'Derry' 'Derry' 'Derry' 'Dory' 'Dory' 'Dory' 'Dory'};
Pen = {'P04' 'P05' 'P08' 'P13' 'P02' 'P03' 'P05' 'P08' 'P00' 'P01' 'P02' 'P04'} ;

for an = 1:length(Animal),



%for g = 1:8
%     Animal = 'Noah';
%     Pen = ['P0' num2str(g)];
    Qualia = 'Good';
%     load(['E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData\' Animal '\Spikes_' Animal '_' Pen '_' Qualia '_Pitch.mat']);
    load(['E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData\' Animal{an} '\Spikes_' Animal{an} '_' Pen{an} '_' Qualia '_Pitch.mat']);

    % stim = 'allHarm';
    % stim2 = 'tone';
    stimList = unique(type);
    Flist = unique(F0);
    UnitList = unique(Y(:,3));
    repList = unique(Y(:,5));
    repList = repList(repList~=0);
    window = [0 0.15]; %seconds

    allResp = [];

    timeCst = diff(window);

    for i = 1:length(UnitList), % for each neural unit
        subY = Y(Y(:,3) == UnitList(i),:);
        allResp(i).Unit = UnitList(i);
        anovaDat = [];
        anovaDat2 = nan(length(Flist)*length(repList),length(stimList));
        gStim = {};
        gF0 = {};
        c = 1;
        for s = 1:length(stimList) % for each stimulus type
            stim = stimList{s};
            for f = 1:length(Flist), % for each F0

                for r = 1:length(repList), %for each repetition/presentation 
                    Nspikes = 0;
                    stimNum = find(strcmp(type,stim) & (F0==Flist(f)));
                    if ~isempty(stimNum),  

                        idx = (subY(:,4) == stimNum) & (subY(:,5)==repList(r)) & (subY(:,2) > window(1)) & (subY(:,2) < window(2));
                        Nspikes = sum(idx);

                        anovaDat(c) = Nspikes ./ timeCst;
                        anovaDat2(f+(r-1)*length(Flist),s) = Nspikes ./ timeCst;
                        gStim{c} = stim;
                        gF0{c} = num2str(Flist(f));
                        c = c+1;
                    end;


                end
            end

        end
        [p,tbl,stats] = anovan(anovaDat',{gStim',gF0'},'model','interaction','display','off','varnames',{'Stimulus' 'F0'});
        allResp(i).p = p;
        allResp(i).tbl = tbl;
        allResp(i).stats = stats;
        allResp(i).values = anovaDat;
        allResp(i).gStim = gStim;
        allResp(i).gF0 = gF0;

        anovaDat3 = anovaDat2(~isnan(sum(anovaDat2,2)),:); % Remove freqs not presented everywhere
        [p,tbl,stats] = anova2(anovaDat3,r,'off');
        allResp(i).anova2.p = p;
        allResp(i).anova2.tbl = tbl;
        allResp(i).anova2.stats = stats;
        allResp(i).anova2.col = stimList;
        allResp(i).anova2.row = Flist(2:end);

    end

%      save(['E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData\' Animal{an} '\Anova2_' Animal{an} '_' Pen{an} '_Pitch.mat'],'allResp','stimList','UnitList','Flist','gStim','gF0');
    save(['E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData\' Animal{an} '\Anova2_' Animal{an} '_' Pen{an} '_Pitch.mat'],'allResp','stimList','UnitList','Flist','gStim','gF0');
end
disp('Done.')

%% One-way ANOVA for Selectivity (for each stim type)

clear;
close all

Animal = {'Ronnie' 'Ronnie' 'Ronnie' 'Ronnie' 'Derry' 'Derry' 'Derry' 'Derry' 'Dory' 'Dory' 'Dory' 'Dory'};
Pen = {'P04' 'P05' 'P08' 'P13' 'P02' 'P03' 'P05' 'P08' 'P00' 'P01' 'P02' 'P04'} ;


% for g = 1:8
for an = 1:length(Animal),
% Animal = 'Noah';
% Pen = ['P0' num2str(g)];
Qualia = 'Good';
% load(['E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData\' Animal '\Spikes_' Animal '_' Pen '_' Qualia '_Pitch.mat']);
load(['E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData\' Animal{an} '\Spikes_' Animal{an} '_' Pen{an} '_' Qualia '_Pitch.mat']);

% load(['E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\Quentin Pitch ephys\PitchEphysData\' 'Spikes_' Animal{an} '_' Pen{an} '_Pitch.mat']);

% stim = 'allHarm';
% stim2 = 'tone';
stimList = unique(type);
Flist = unique(F0);
UnitList = unique(Y(:,3));
repList = unique(Y(:,5));
repList = repList(repList~=0);
window = [0 0.15]; %seconds

allResp = [];

timeCst = diff(window);

for i = 1:length(UnitList),
    subY = Y(Y(:,3) == UnitList(i),:);
    allResp(i).Unit = UnitList(i);
    
    for s = 1:length(stimList),
        stim = stimList{s};
        anovaDat = [];
        gF0 = {};
        c = 1;
        for f = 1:length(Flist),
            
            for r = 1:length(repList),
                Nspikes = 0;
                stimNum = find(strcmp(type,stim) & (F0==Flist(f)));
                if ~isempty(stimNum),  
                
                    idx = (subY(:,4) == stimNum) & (subY(:,5)==repList(r)) & (subY(:,2) > window(1)) & (subY(:,2) < window(2));
                    Nspikes = sum(idx);
                    anovaDat(c) = Nspikes ./ timeCst;
                    gF0{c} = num2str(Flist(f));
                    c = c+1;
                end;
                
            end
        end
        
        [p(s),tbl{s},stats{s}] = anova1(anovaDat,gF0,'off'); %do a 1-way anova for effects of F0 on spike rate for this stimulus type
        anovaDat2{s} = anovaDat;
        
    end
    allResp(i).p = p; %one-way anova on this unit (effect of F0 on spike rate), for each stim type
    allResp(i).tbl = tbl;
    allResp(i).stats = stats;
    allResp(i).values = anovaDat2; %the spike rates for each trial that we put into the 1-way anova
    
end

save(['Anova1_' Animal{an} '_' Pen{an} '_Pitch.mat'],'allResp','stimList','UnitList','Flist','gF0');
% save(['Anova1_' Animal '_' Pen '_Pitch.mat'],'allResp','stimList','UnitList','Flist','gF0');
end
disp('Done.')

%% Unit Selection
%what is this section? Can it be run then forgotten about?
clear;
close all;
cd('E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData')

%fList = {'a1_P00_Pitch' 'a2_P00_Pitch' 'Ronnie_P04_Pitch' 'Ronnie_P13_Pitch' 'Derry_P02_Pitch' 'Derry_P03_Pitch'};
fList = {'Derry_P02_Pitch' 'Derry_P03_Pitch' 'Derry_P05_Pitch' 'Derry_P08_Pitch' 'Ronnie_P04_Pitch' 'Ronnie_P05_Pitch' 'Ronnie_P08_Pitch' 'Ronnie_P13_Pitch' 'Dory_P00_Pitch' 'Dory_P01_Pitch' 'Dory_P02_Pitch' 'Dory_P04_Pitch'...
'Noah_P01_Pitch' 'Noah_P02_Pitch' 'Noah_P03_Pitch' 'Noah_P04_Pitch' 'Noah_P05_Pitch' 'Noah_P06_Pitch' 'Noah_P07_Pitch' 'Noah_P08_Pitch'};


allDat = {}; %blank cell array

for i = 1:length(fList), %for each penetration in each animal 
    fName = fList{i}; %creates a variable for whichever file is to be loaded this time
    
    load(['E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData\Responsiveness_' fName '.mat'])
    allDat{i}.ttestPval = allResp; %ttest on response of neuron during sound vs silence (at each 
    
    
    load(['E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData\Anova1_' fName '.mat'])
    allDat{i}.DataAnova1 = allResp; %effect of F0 on spike rate, for each of the stimuli types
    
    load(['E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData\Anova2_' fName '.mat'])
    allDat{i}.Data = allResp;
    allDat{i}.Units = UnitList;
    allDat{i}.Animal(1:length(allResp)) = 1;
    allDat{i}.fName = fName;
    
    
    % Select responsive units (Anova and interaction < 0.05)
    allDat{i}.idxRespStim = arrayfun(@(x)(x.p(1)<0.05),allDat{i}.Data);
    allDat{i}.idxRespF0 = arrayfun(@(x)(x.p(2)<0.05),allDat{i}.Data);
    allDat{i}.idxRespInteraction = arrayfun(@(x)(x.p(3)<0.05),allDat{i}.Data);
    
    % ttest
    allDat{i}.idxRespTtest = arrayfun(@(x)(sum(x.Global.p<0.01)>0),allDat{i}.ttestPval);
    
    % Anova1 (pass test for all Harm, p<0.05)
    allDat{i}.idxRespAnova1 = arrayfun(@(x)(x.p(8)<0.05),allDat{i}.DataAnova1);
    
    % Anova1 hard (pass test for all Harm, p<0.0038 - Bonferronni
    % correctionn)
    allDat{i}.idxRespAnova1BC = arrayfun(@(x)(x.p(8)<0.0038),allDat{i}.DataAnova1);
    
    
%     allDat{i}.idxResp = (allDat{i}.idxRespF0 | allDat{i}.idxRespInteraction) & allDat{i}.idxRespTtest;
 %   allDat{i}.idxResp = allDat{i}.idxRespAnova1 & allDat{i}.idxRespTtest;
    allDat{i}.idxResp = allDat{i}.idxRespF0 | allDat{i}.idxRespInteraction; %Did the 2-way anova show a main effect of F0 or a F0/stimtype interaction?

end

save(['allDat_allanimalsandpens.mat'],'allDat','stimList', '-v7.3');
disp('Done.')
% [Animal(idxResp)' Units(idxResp)];
% 
% [c,m] = multcompare(Data(1).stats,'Alpha',0.05,'CType','tukey-kramer','Display','off','Dimension',[1]);
% idxSigDiff = ~(c(:,3)<0 & c(:,5)>0);
% c(idxSigDiff,[1 2]);

%% The main population analyses for the SfN 2019 poster are below %%

%% Totals
% total number of units, and number responsive and tuned (2-way Anova):
% clear;
% close all;

cd('E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData');
% load('allDat_allanimalsandpens.mat'); %load all the data with Noah

% cd('E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys');
% load('allDat_allanimalsandpens_noNoah.mat'); %load all the data without Noah (sanity check)

numunits=0; %set number of units to 0
p_ttest=[]; %create an array for the significance of t test results
tuned=[]; %create an array for the number of tuned neurons
for ii=1:length(allDat) %for all the data
    numunits=numunits+(length(allDat{ii}.ttestPval)); %add unit to total number if it is has a t-test value (should be every neuron)
    d=allDat{ii}.ttestPval; %this is every neuron with a t-test value
    for jj=1:length(d) %for every neuron 
        p_ttest=[p_ttest; d(1).Global.p]; %results of t-test for said neuron 
    end
    tuned=[tuned allDat{ii}.idxResp]; %number of neurons that are sound-responsive
end
disp(numunits) %total number of units 
disp(sum(p_ttest<0.05)) %number of units that had a signif different spike rate (paired t-test) during sound and silence (pooled across all F0 and stim types)
disp(sum(tuned)) %number of tuned neurons

%% Plot tuning curve for the best tuned units
% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13
% %load(['/Users/bras1912/Google Drive/Quentin Pitch ephys/PitchEphysData/Anova2_Derry_P02_Pitch.mat']); %just to get stimList
% clear;
% close all;
% 
%load(['E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\Quentin Pitch ephys\PitchEphysData\Anova2_Derry_P02_Pitch.mat']); %just to get stimList
cd ('E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData')
load(['allDat_allanimalsandpens.mat']); %load the master data file
%Make a list of all 1-way ANOVA p-values, and
%the corresponding file indices
allp=[]; allpen=[]; allunits=[]; %sets up arrays for all total p values, penetrations, and unit numbers
for ii=1:length(allDat) %for each file (penetration)
    for jj=1:length(allDat{ii}.DataAnova1) %for each unit that is responsive to F0          
            allp=[allp; allDat{ii}.DataAnova1(jj).p;]; %p values for all 13 stimuli %during this loop, keep adding new p values to this array)
            allpen=[allpen; ii;]; %penetration index         
            allunits=[allunits; jj;]; %units index         
    end
end

save('All1wayanova.mat','allDat','allp','allpen','allunits', '-v7.3'); %save what has been done

stim='tone'; %your stimulus of interest, from the list in the comment above
stimind=find(strcmp(stimList,stim)); %your stimulus of interest, from the list in the comment above
fsig=find(allp(:,stimind)<0.01); %finds all the significant units for that stimulus
% 
% %[a,b]=sort(allp(:,stimind));
% %hist(a,100)
% % %reorder according to p-value of anova
% % allp=allp(b,:);
% % allpen=allpen(b);
% % allunits=allunits(b);
% 
% disp("Done")


for ii=1:length(fsig) %for each significant unit
%for ii=1:5 %for each significant unit
    thisunit=fsig(ii); %create a variable
    spikerate=allDat{allpen(thisunit)}.Data(allunits(thisunit)).values; %spike rate value
    allstims=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gStim; %stimulus value
    f0=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gF0; %frequency value
    f0=cellfun(@str2num, f0); %make F0 values numeric
    stim2use=strcmp(allstims,stim); %find trials with this stimulus type
    
    trials=find(stim2use);
    f0=f0(trials); %frequency is equal to the appropriate one
    allf0=unique(f0); %finds all the different frequencies for above
    spikerate=spikerate(trials); %finds the spike rate for trials
    mrate=[]; %array ready for mean spike rate
    serate=[]; %array ready for standard error of spike rate
    for jj=1:length(allf0) %for every frequency
        mrate(jj)=mean(spikerate(find(f0==allf0(jj)))); %calculate mean spike rate
        serate(jj)=std(spikerate(find(f0==allf0(jj))))/sqrt(length(find(f0==allf0(jj)))); %calculate standard error for spike rate
    end
    figure(1); %make a figure
    clf %clear previous stuff
    shadedErrorBar(log(allf0),mrate,serate,'-b'); 
    set(gca,'xtick',log(allf0(1:3:end)),'xticklabels',(allf0(1:3:end)),'fontname','ariel','fontsize',16);
    xlabel('F0 (Hz)','fontname','ariel','fontsize',20);
    ylabel('Spike rate (Hz)','fontname','ariel','fontsize',20);
    xlim([log(allf0(1)) log(allf0(end))]);
    title(['Stim = ' stim '; Unit = ' num2str(thisunit) '; p = ' num2str(allp(thisunit,stimind))]);
    box off
    ylim([0 40])
    set(gcf,'position',[500 466 477 339])
    set(gcf,'paperpositionmode','auto')
    pstr=num2str(allp(thisunit,stimind));
    
%     print([pwd '/Figures/' 'Tuning/' stim '_' pstr(3:end) '_' num2str(thisunit) '.png'],'-dpng')
    

%     pause()
    
end
 
 %% for plotting tuning curve of one unit of your choice
 %first run the above section
%stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
%             1       2          3         4        5             6          7                 8           9          10       11       12        13               
%For Noah: 'CTO'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'SAMtones'  'allHarm'  'allHarmRand' 'alt'  'high'    'low'   'rand'  'tone'      

% load('E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\Quentin Pitch ephys\PitchEphysData\Spikes_Ronnie_P05_Pitch.mat');
Animal = 'Noah';
Pen = 'P01';
Qualia = 'Good';
load(['E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData\' Animal '\Spikes_' Animal '_' Pen '_' Qualia '_Pitch.mat']);
stimList = unique(type);
stim='CT0'; %your stimulus of interest, from the list in the comment above
stimind=find(strcmp(stimList,stim)); %your stimulus of interest, from the list in the comment above
f1 = figure;

for ii=105 %pick a unit between 1 and 1143
    thisunit=ii; %used to get data for this unit
    spikerate=allDat{allpen(thisunit)}.Data(allunits(thisunit)).values;
    allstims=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gStim;
    f0=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gF0;
    f0=cellfun(@str2num, f0); %make F0 values numeric
    stim2use=strcmp(allstims,stim); %find trials with this stimulus type
    
    trials=find(stim2use); %get data for the correct trials
    f0=f0(trials); %possible frequencies
    allf0=unique(f0); %list all the different frequencies used 
    spikerate=spikerate(trials); %spike rates for this trial
    mrate=[]; 
    serate=[];
    for jj=1:length(allf0) %for each frequency tested
        mrate(jj)=mean(spikerate(find(f0==allf0(jj)))); %what was the mean spike rate
        serate(jj)=std(spikerate(find(f0==allf0(jj))))/sqrt(length(find(f0==allf0(jj)))); %what was the SE of the spike rate
    end
    figure(f1);
    clf
    shadedErrorBar(log(allf0),mrate,serate,'-b');
    set(gca,'xtick',log(allf0(1:4:end)),'xticklabels',(allf0(1:4:end)),'fontname','ariel','fontsize',16);
    set(gcf, 'DefaultFigureVisible', 'on')
    xlabel('F0 (Hz)','fontname','ariel','fontsize',26);
    ylabel('Spike rate (Hz)','fontname','ariel','fontsize',26);
    xlim([log(allf0(1)) log(allf0(end))]);
    %title(['Stim = ' stim '; Unit = ' num2str(thisunit) '; p = ' num2str(allp(thisunit,stimind))]);
     box off
    ylim([0 50])
    set(gcf,'position',[500 466 477 339])
    set(gcf,'paperpositionmode','auto')
    pstr=num2str(allp(thisunit));
    
    %print([pwd '/Figures/' 'Tuning/' stim '_' pstr(3:end) '_' num2str(thisunit) '.png'],'-dpng')
    

    
    
end

%% Plot tuning curves for multiple stimuli 
%Run unit selection and best tuned units first
%stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
%             1       2          3         4        5             6          7                 8           9          10       11       12        13
%plotstims={'CT0'; 'low'; 'F0MaskLow'; 'high'; F0MaskLow}; %your stimulus of interest, from the list in the comment above
%cc={'k','b','c','r','m'}; %colors
%plotstims={'CT0'; 'F0MaskLow'; 'F0MaskHigh'}; %your stimulus of interest, from the list in the comment above
%cc={'k','c','m'}; %colors
% plotstims={'CT0'; 'low'; 'high'}; %your stimulus of interest, from the list in the comment above
% cc={'k','b','r'}; %colors
% plotstims={'CT0'; 'CT5'; 'CT10'; 'CT20';'CT40'};
% cc={'k','b','c','m','r'}; %colors
plotstims={'CT0'; 'tone'};
cc={'k','g','r'}; %colors
% plotstims={'CT0'; 'low'; 'high'}; %tuning curves of these 2 stimuli
% cc={'k','b','r'}; %colors

stim1='CT0'; %your stimulus of interest, from the list in the comment above
stimind1=find(strcmp(stimList,stim1)); %your stimulus of interest, from the list in the comment above
fsig1=find(allp(:,stimind1)<0.05); 
stim2='tone';
stimind=find(strcmp(stimList,stim2)); %your stimulus of interest, from the list in the comment above
fsig2=find(allp(:,stimind)>0.05);
% stim3='rand'; %third stimulus 
% stimind3=find(strcmp(stimList,stim3));
% fsig3=find(allp(:,stimind)>0.05);

% fsig=intersect(fsig1,fsig2);
% fsigall=intersect(fsig,fsig3);
% fsig=fsig1;


for ii=1:length(fsig1), %for each significant unit
    figure(1);
    clf
    hold on
    thisunit=fsig1(ii);
    spikerate=allDat{allpen(thisunit)}.Data(allunits(thisunit)).values;
    allstims=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gStim;
    f0=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gF0;
    f0=cellfun(@str2num, f0); %make F0 values numeric
    
    bestf0=[];
    for kk=1:length(plotstims)
        spikerate=allDat{allpen(thisunit)}.Data(allunits(thisunit)).values;
        allstims=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gStim;
        f0=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gF0;
        f0=cellfun(@str2num, f0); %make F0 values numeric
        
        stim=plotstims{kk};
        stim2use=strcmp(allstims,stim); %find trials with this stimulus type

        trials=find(stim2use);
        f0=f0(trials);
        allf0=unique(f0);
        spikerate=spikerate(trials);
        mrate=[];
        serate=[];
        for jj=1:length(allf0)
            mrate(jj)=mean(spikerate(find(f0==allf0(jj))));
            serate(jj)=std(spikerate(find(f0==allf0(jj))))/sqrt(length(find(f0==allf0(jj))));
        end
        figure(1);
        shadedErrorBar(log(allf0),mrate,serate,cc{kk},1);
        %errorbar(log(allf0),mrate,serate,cc{kk});
        %plot(log(allf0),mrate,cc{kk});
        [temp,bestf0ind]=max(mrate);
        bestf0(kk)=allf0(bestf0ind); %save the peak of the tuning curve here
    end
    set(gca,'xtick',log(allf0(1:4:end)),'xticklabels',(allf0(1:4:end)),'fontname','ariel','fontsize',20);
    xlabel('F0 (Hz)','fontname','ariel','fontsize',22);
    ylabel('Spike rate (Hz)','fontname','ariel','fontsize',22);
    xlim([log(allf0(1)) log(allf0(end))]);
%     title(num2str(bestf0));
%     title([stim1,' vs ',stim2]);
    box off
    %ylim([0 40])
    set(gcf,'position',[794   459   477   339])
    set(gcf,'paperpositionmode','auto')
    pstr=num2str(allp(thisunit,stimind));
    
%     legend('','','',stim1,'','','',stim2,'','','',stim3);
    
%     print([pwd '/Figures/' 'Tuning/MultTuning_' stim1 '_' stim2 '_' stim3 '_' num2str(thisunit) '.png'],'-dpng')    
    print([pwd '/Figures/' 'Tuning/MultTuning_' stim1 '_' stim2 '_' num2str(thisunit) '.png'],'-dpng')
%      print([pwd '/Figures/' 'Tuning/MultTuning_' stim1 '_'  num2str(thisunit) '.png'],'-dpng')


%      pause()
    
end

%% Plot % neurons FO-sensitive for different stim types
%Run unit selection and best tuning curves first
% cd ('E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData')
% load('All1wayanova.mat');
% stims={'CT0'; 'tone'; 'low'; 'F0MaskLow'; 'high'; 'F0MaskHigh'}; 
% stims={'CT0'; 'CT5'; 'CT10'; 'CT20';'CT40'}; 
stims={'high'; 'alt'};
percenttuned=[];
for ii=1:length(stims) %for every stimulus
    stimind=find(strcmp(stimList,stims{ii})); %find the stimulus in the list
    fsig=find(allp(:,stimind)<0.05); %find all the units tuned to thist stimulus
    percenttuned(ii)=length(fsig)/size(allp,1)*100; %find the percentage tuned from these
end

% lowstim=find(strcmp(stimList,stims{3}));
% fsiglow=find(allp(:,lowstim)<0.05);
% highstim=find(strcmp(stimList,stims{5}));
% fsighigh=find(allp(:,highstim)<0.05);
% doublesig=intersect(fsiglow,fsighigh);
% percentlh=length(doublesig)/size(allp,1)*100;
% percenttuned(7)=percentlh;

figure(2);clf;hold on; %make a figure window
bar([1:length(stims)]*4,percenttuned) %make a bar chart for all stim on x, against % tuned
set(gca,'xtick',[1:length(stims)]*4,'xticklabels',stims,'fontname','ariel','fontsize',40);
%xlabel('F0 (Hz)','fontname','ariel','fontsize',32);
ylabel('F0-sensitive neurons (%)','fontname','ariel','fontsize',40);
box off

%now find percent also tuned to a control stimulus:
% stims2={''; 'tone'; 'low'; 'F0MaskLow'; 'high'; 'F0MaskHigh'}; %your stimulus of interest, from the list in the comment above
stims2={''; 'CT5'; 'CT10'; 'CT20';'CT40'}; 
controlstim='CT0'; 
stimind2=find(strcmp(stimList,controlstim)); 
fsig_control=find(allp(:,stimind2)<0.05); %indices of neurons tuned to F0 for the control stim
percenttuned_both=[];
for ii=1:length(stims2)
    stimind3=find(strcmp(stimList,stims2{ii})); %your stimulus of interest, from the list in the comment above
    fsig3=find(allp(:,stimind3)<0.05); %index of all neurons tuned to F0 for other stims
    fsigboth=intersect(fsig3,fsig_control); %all units tuned to both
    percenttuned_both(ii)=length(fsigboth)/size(allp,1)*100; 
end

% fsigall=intersect(doublesig,fsig_control);
% percenttuned_both(7)=length(fsigall)/size(allp,1)*100;

% figure(3);clf;hold on;
bar([1:length(stims2)]*4,percenttuned_both,'r')
xlim([2 length(stims2)*4+6])
% ylabel({'% neurons F0-sensitive';'for unjittered clicks'},'fontname','ariel','fontsize',20);
% print([pwd '/Figures/' 'Tuning/PercentTuned.png'],'-dpng')

hold off
%set(gcf,'position',[582   465   364   317])
set(gcf,'position',[582   503   781   279])
set(gcf,'paperpositionmode','manual')
set(gcf,'paperposition',[0 0 30 25])
%set(gca,'xtick',[1:length(stims)+1]*4,'xticklabels',{'0';'5';'10';'20';'40'},'fontname','ariel','fontsize',20);
%xlabel('% jitter')
% print([pwd '/Figures/' 'Tuning/PercentTuned_jitter.png'],'-dpng')
print([pwd '/Figures/' 'Tuning/PercentTuned_talk.png'],'-dpng')

%% Plot responses to one stimuli against another
% Run Unit Selection first
close all
fList = {'Derry_P02_Pitch' 'Derry_P03_Pitch' 'Derry_P05_Pitch' 'Derry_P08_Pitch'...
    'Ronnie_P04_Pitch' 'Ronnie_P05_Pitch' 'Ronnie_P08_Pitch' 'Ronnie_P13_Pitch'...
    'Dory_P00_Pitch' 'Dory_P01_Pitch' 'Dory_P02_Pitch' 'Dory_P04_Pitch'...
    'Noah_P01_Pitch' 'Noah_P02_Pitch' 'Noah_P03_Pitch' 'Noah_P04_Pitch' ...
    'Noah_P05_Pitch' 'Noah_P06_Pitch' 'Noah_P07_Pitch' 'Noah_P08_Pitch' };
CFValues =[];
CFidx=[];
%'CT0' 'CT10' 'CT20' 'CT40' 'CT5' 'F0MaskHigh' 'F0MaskLow' 'allHarm' 'alt'  'high'  'low'  'rand'  'tone'
allstims={'CT0';'CT10'; 'CT20'; 'CT40'; 'CT5'; 'F0MaskHigh'; 'F0MaskLow'; 'allHarm'; 'alt';  'high'; 'low'; 'rand';  'tone'};
for ii = 1:length(fList), %for each penetration
    uIdx = find(allDat{ii}.idxResp);
    %CFValues = zeros(length(uIdx),2);

    % load anova2 to get the right gstim
    fName = fList{ii};
    load(['E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData\Anova2_' fName '.mat'])
    stim1='rand'; %y axis
    stim2='high'; %x axis
    stim1Idx = strcmp(gStim,stim1);
    stim2Idx = strcmp(gStim,stim2);
    

    for i = 1:length(uIdx), %for each neuron
        if sum(allResp(uIdx(i)).values(stim1Idx))>0 && sum(allResp(uIdx(i)).values(stim2Idx))>0 %if some spikes for both stimuli
            val = allResp(uIdx(i)).values(stim1Idx);
            Freq = str2double(gF0(stim1Idx));
            FList = unique(Freq);
            meanMat=[];
            errMat=[];
            for j = 1:length(FList),
                meanMat(j) = mean(val(FList(j)==Freq));
                errMat(j) = ste(val(FList(j)==Freq));
            end
            [~, iV] = max(meanMat); %get the largest mean value for stim1 (which frequency has best resp)
            CFV1=FList(iV);
            

            val = allResp(uIdx(i)).values(stim2Idx);
            Freq = str2double(gF0(stim2Idx));
            FList = unique(Freq);
            meanMat=[];
            errMat=[];
            for j = 1:length(FList),
                meanMat(j) = mean(val(FList(j)==Freq));
                errMat(j) = ste(val(FList(j)==Freq));
            end
            [~, iV] = max(meanMat); %get the largest mean value for stim2
            CFV2=FList(iV);
            
            CFValues = [CFValues; CFV1 CFV2]; %mean value array for each unit
        end

    end
end

CFidx(:,1) = arrayfun(@(x)(find(x==FList)),CFValues(:,1)); %stim 1 CF index
CFidx(:,2) = arrayfun(@(x)(find(x==FList)),CFValues(:,2)); %stim 2 CF index
C=confusionmat(CFidx(:,1),CFidx(:,2));


if ((sum(C(1,:))==0) || (sum(C(:,1))==0)), %this happens when low harm is used, as lowest 2 F0s are missing
    Csm=C(3:end,3:end);
    MI = MI2fromtable(Csm) %Mutual Information
    disp('^^Mutual information')
    QOctave = (sum(diag(Csm))+sum(diag(Csm,-1))+sum(diag(Csm,1)))/(sum(sum(Csm)))*100; % percent accuracy within 1/4 octave
    disp('^^%within 1/4 octave')
else
   MI = MI2fromtable(C) %Mutual Information
   disp('^^Mutual information')
   QOctave = (sum(diag(C))+sum(diag(C,-1))+sum(diag(C,1)))/(sum(sum(C)))*100 % percent accuracy within 1/4 octave
   disp('^^%within 1/4 octave')
end



% figure(1);
% scatterQ(CFValues(:,1),CFValues(:,2));
% [r, p] = corrcoef(CFValues(:,1),CFValues(:,2));
% title(['r = ' num2str(r(2,1)) ', p = ' num2str(p(2,1))]);
% 
% CFidx(:,1) = arrayfun(@(x)(find(x==FList)),CFValues(:,1));
% CFidx(:,2) = arrayfun(@(x)(find(x==FList)),CFValues(:,2));

% figure(2);
% idxDiff = CFidx(:,1) - CFidx(:,2);
% hist(idxDiff,25)

% figure(3);
% octDiff = log2(CFValues(:,1)./CFValues(:,2));
% hist(octDiff,25)
% xlabel('CF distance (octaves)')

% figure(4);
% fDiff = CFValues(:,1) - CFValues(:,2);
% hist(fDiff,20)

figure(5);

C=confusionmat(CFidx(:,1),CFidx(:,2));
imagesc(flipud(C));
axis square
freqs=unique(FList);
set(gca,'ytick',[1:4:length(FList)],'yticklabels',fliplr(freqs(1:4:end)),'fontname','ariel','fontsize',42)
set(gca,'xtick',[1:4:length(FList)],'xticklabels',freqs(1:4:end))
xlabel(stim2,'fontsize',52)
ylabel(stim1,'fontsize',52)
colorbar
set(gcf,'position',[803   362   560   420])
set(gcf,'paperpositionmode','manual')
set(gcf,'paperposition', [0 0 40 25])
% title({['Mutual Information = ' num2str(MI)],['% within 1/4 octave =' num2str(QOctave)]});
% ax.FontSize = 12;
% ax.TitleFontSizeMultiplier = 0.5;
print([pwd '/Figures/' 'Tuning/confusionTuning_' stim1 '_' stim2 '.png'],'-dpng','-r0')

   
  %  pause()
  
 (sum(sum(C)))

%% Find putative pitch neurons
% clear all
cd('E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData')
% load('All1wayanova.mat')
% load('allDat_allanimalsandpens.mat')
%stims={'CT0';'CT10'; 'CT20'; 'CT40'; 'CT5'; 'F0MaskHigh'; 'F0MaskLow'; 'allHarm'; 'alt';  'high'; 'low'; 'rand';  'tone'};
stims_tuned={'low';'F0MaskLow';'high'; 'F0MaskHigh';'tone';'CT0';'CT5';'CT10'}; %want tuning for these
stims_untuned={'CT20';'CT40';'rand';'alt'}; %original - don't want tuning

%No matter what you put in the untuned box, it only cares about the second
%one... -Dan

cc_t={'b','c','r','m','g','k','k','k','y','y'};
cc_u={'k','k','r','b'};
units_tuned=zeros(size(allp,1),length(stims_tuned)); %initialise tuned array
units_untuned=zeros(size(allp,1),length(stims_untuned)); %initialise untuned array


for ii=1:length(stims_tuned)
    stimind=find(strcmp(stimList,stims_tuned{ii})); %your stimulus of interest, from the list in the comment above
    units_tuned(find(allp(:,stimind)<0.01),ii)=1; %all units tuned for this stimuli
end
for ii=1:length(stims_untuned)
    stimind=find(strcmp(stimList,stims_untuned{ii})); %your stimulus of interest, from the list in the comment above
    units_untuned(find(allp(:,stimind)>0.01),ii)=1; %all units untuned for this stimuli
    %units_untuned(ii).units=find(allp(:,stimind)>0.05);
end

%pitch_units=intersect(find(sum(units_tuned,2)>(length(stims_tuned)-2)),find(sum(units_untuned,2)>(length(stims_untuned)-2)));
%pitch_units=find(sum(units_tuned,2)>2);
pitch_units=intersect(find(sum(units_tuned,2)>4),find(units_untuned(:,2)>0)); %original - at least 3 pitch stimuli but not CT40
% pitch_units=intersect(find(sum(units_tuned,2)>4),find(sum(units_untuned,2)>1)); %find neurons that have tuning for at least 5 pitch stimuli but and at least one untuned

    
length(pitch_units) %how many pitch neurons are there

for ii=1:length(pitch_units) %for every putative neuron
    allmrate=[];
    thisunit=pitch_units(ii);
    figure(1);hold on;clf
    bestf0=[];
    for kk=1:length(stims_tuned) %for every tuned stimulus
        spikerate=allDat{allpen(thisunit)}.Data(allunits(thisunit)).values; %find spike rate
        allstims=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gStim; %find stimulus
        f0=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gF0; %find frequencies
        f0=cellfun(@str2num, f0); %make F0 values numeric
        
        stim=stims_tuned{kk};
        stim2use=strcmp(allstims,stim); %find trials with this stimulus type

        trials=find(stim2use);
        f0=f0(trials);
        allf0=unique(f0); %different frequencies of all valid trials
        spikerate=spikerate(trials); %spike rate of valid trials
        mrate=[];
        serate=[];
        for jj=1:length(allf0) %for every frequency
            mrate(jj)=mean(spikerate(find(f0==allf0(jj)))); %calculate the mean spike rate
            serate(jj)=std(spikerate(find(f0==allf0(jj))))/sqrt(length(find(f0==allf0(jj)))); %calculate the SE
        end
        if strcmp(stim,'low') %if the stimulus type is low harmonics
            allmrate=[allmrate; [0 0 mrate]];
        elseif strcmp(stim,'F0MaskLow')
            allmrate=[allmrate; [0 0 mrate]];
        else
            allmrate=[allmrate; mrate];
        end
        figure(1);hold on
        
        %shadedErrorBar(log(allf0),mrate,serate,cc{kk},1);
        %errorbar(log(allf0),mrate,serate,cc{kk});
        
        if strcmp(stim,'CT5')
            plot(log(allf0),mrate,'color',[.5 .5 .5],'linestyle','-');
        else
            plot(log(allf0),mrate,cc_t{kk},'linestyle','-');
        end
        [temp,bestf0ind]=max(mrate);
        bestf0(kk)=allf0(bestf0ind); %save the peak of the tuning curve here
    end
    for kk=1:length(stims_untuned)
        spikerate=allDat{allpen(thisunit)}.Data(allunits(thisunit)).values;
        allstims=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gStim;
        f0=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gF0;
        f0=cellfun(@str2num, f0); %make F0 values numeric
        
        stim=stims_untuned{kk};
        stim2use=strcmp(allstims,stim); %find trials with this stimulus type

        trials=find(stim2use);
        f0=f0(trials);
        allf0=unique(f0);
        spikerate=spikerate(trials);
        mrate=[];
        serate=[];
        for jj=1:length(allf0)
            mrate(jj)=mean(spikerate(find(f0==allf0(jj))));
            serate(jj)=std(spikerate(find(f0==allf0(jj))))/sqrt(length(find(f0==allf0(jj))));
        end
        allmrate=[allmrate; mrate];
        figure(1);hold on
        %shadedErrorBar(log(allf0),mrate,serate,cc{kk},1);
        %errorbar(log(allf0),mrate,serate,cc{kk});
        plot(log(allf0),mrate,cc_u{kk},'linestyle','--');
    end
    set(gca,'xtick',log(allf0(1:3:end)),'xticklabels',(allf0(1:3:end)),'fontname','ariel','fontsize',24);
    xlabel('F0 (Hz)','fontname','ariel','fontsize',28);
    ylabel('Spike rate (Hz)','fontname','ariel','fontsize',28);
    xlim([log(allf0(1)) log(allf0(end))]);
    title(num2str(bestf0));
    box off
    %ylim([0 40])
    set(gcf,'position',[794   459   477   339])
    set(gcf,'paperpositionmode','auto')
    pstr=num2str(allp(thisunit,stimind));   
    print([pwd '/Figures/' 'Tuning/PitchNeuron_' num2str(thisunit) '.png'],'-dpng')

    figure(2);clf
    imagesc(allmrate);
    set(gca,'xtick',[1:4:length(allf0)],'xticklabels',(allf0(1:4:end)),'fontname','ariel','fontsize',40);
    stimstuneduntuned=[stims_tuned;stims_untuned];
    set(gca,'ytick',[1:length(stimstuneduntuned)],'yticklabels',stimstuneduntuned)
    colorbar
%     title(allDat{allpen(thisunit)}.fName)
    set(gcf,'position',[679   459   592   339])
    set(gcf,'paperpositionmode','manual')
    set(gcf,'paperposition', [0 0 30 20])
    print([pwd '/Figures/' 'Tuning/PitchNeuron_' num2str(thisunit) '_imagesc.png'],'-dpng')
    
end

save('pitch_information', 'allmrate', 'stimstuneduntuned', 'pitch_units')

%% Pitch neuron locations 


% fList = {'Derry_P02_Pitch' 'Derry_P03_Pitch' 'Derry_P05_Pitch' 'Derry_P08_Pitch' 'Ronnie_P04_Pitch' 'Ronnie_P05_Pitch' 'Ronnie_P08_Pitch' 'Ronnie_P13_Pitch' 'Dory_P00_Pitch' 'Dory_P01_Pitch' 'Dory_P02_Pitch' 'Dory_P04_Pitch'};

fList = {'Derry_P02_Pitch' 'Derry_P03_Pitch' 'Derry_P05_Pitch' 'Derry_P08_Pitch'...
    'Ronnie_P04_Pitch' 'Ronnie_P05_Pitch' 'Ronnie_P08_Pitch' 'Ronnie_P13_Pitch'...
    'Dory_P00_Pitch' 'Dory_P01_Pitch' 'Dory_P02_Pitch' 'Dory_P04_Pitch'...
    'Noah_P01_Pitch' 'Noah_P02_Pitch' 'Noah_P03_Pitch' 'Noah_P04_Pitch' ...
    'Noah_P05_Pitch' 'Noah_P06_Pitch' 'Noah_P07_Pitch' 'Noah_P08_Pitch' };

%fList = {'Ronnie_P05_Pitch' 'Ronnie_P08_Pitch' 'Ronnie_P04_Pitch'...
 %   'Ronnie_P13_Pitch' 'Derry_P02_Pitch' 'Derry_P03_Pitch' ...
 %   'Derry_P05_Pitch' 'Derry_P08_Pitch' ...
 %   'Dory_P00_Pitch' 'Dory_P01_Pitch' 'Dory_P02_Pitch' 'Dory_P04_Pitch'};
%LocPen = {'LowA1' 'LowAAF' 'LowAAF' 'MidA1' 'LowAAF' 'MidA1' 'LowA1' 'LowA1' ...
 %   'LowA1' 'LowAAF' 'MidAAF' 'LowAAF'};
LocPen = {'LowAAF' 'MidA1' 'LowA1' 'LowA1'...
    'LowAAF' 'LowAAF' 'LowA1' 'MidA1'...
    'LowA1' 'LowAAF' 'MidAAF' 'LowAAF'...
    ' ' 'LowA1' 'HighA1' 'LowA1' ...
    'LowPPF' 'LowAAF' 'LowPPF' 'High AAF'};

LocPenidx = [2 3 1 1 2 2 1 3 1 2 4 2 1 2 5 1 7 2 7 6];

%1 = LowA1, 2 = LowAAF, 3 = MidA1, 4 = MidAAF, 5 = HighA1, 6 = High AAF,
% 7 = LowPPF

PListidx = unique(LocPenidx);
PList = unique(LocPen);
finalDat = [];
clear loc_pitchUnits fname_pitchUnits
for pp=1:length(pitch_units)
    thisunit=pitch_units(pp);
    Dat = allDat{allpen(thisunit)};
    loc_pitchUnits{pp}=LocPen{allpen(thisunit)};
    fname_pitchUnits{pp}=fList{allpen(thisunit)};
end
  
fname_pitchUnits
loc_pitchUnits 

Pcount_all=zeros(1,length(PList));
Pcount_pitch=zeros(1,length(PList));
for pp=1:length(allDat)
    numu=length(allDat{pp}.ttestPval);%num units in this pen
    Pcount_all(LocPenidx(pp))=Pcount_all(LocPenidx(pp))+numu;
end
for pp=1:length(pitch_units)
    thisunit=pitch_units(pp);
    thisarea = LocPenidx(allpen(thisunit));
    Pcount_pitch(thisarea)=Pcount_pitch(thisarea)+1;
end
 
PList
Pcount_all
Pcount_pitch
Pcount_pitch./Pcount_all*100
% for imp = 1:length(allDat),
%     Dat = allDat{imp};
%     uIdx = find(Dat.idxResp);
%     CFValues = zeros(length(uIdx),2);
%     nU = length(uIdx);
%     nAllU = length(Dat.idxResp);
%     nUTones = sum(arrayfun(@(x)(x.p(13)<0.05),Dat.DataAnova1));
%     nUResp = sum(Dat.idxRespTtest);
%     nUBoth = sum(arrayfun(@(x)(x.p(13)<0.05),Dat.DataAnova1) & Dat.idxResp);
%     
%     finalDat(imp,:) = [nAllU nUResp nU nUTones nUBoth];
%     
% end
% finalDat2 = [];
% finalDat2prc = [];
% for i = 1:length(PList)
%     idx = strcmp(LocPen,PList(i));
%     finalDat2(i,:) = sum(finalDat(idx,:),1);
%     finalDat2prc(i,:) = round(finalDat2(i,:) ./ finalDat2(i,1) *100,2);
% end
% disp(finalDat2)
% disp(finalDat2prc)
% disp(PList)
%% Find temporal encoding neurons
cd('E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData')
% load('All1wayanova.mat')
% load('allDat_allanimalsandpens.mat')

stims_tuned={'high';'alt';'CT0';'F0MaskHigh';'tone'}; %want tuning for these, make sure alt and rand are the last 2 and the 2 before them are normal!
stims_untuned={'CT40'}; %don't want tuning for these
stimsforaxes={'low';'F0MaskLow';'high';'F0MaskHigh';'tone';'CT0';'CT5';'CT10';'CT20';'CT40';'rand';'alt'};

cc_t={'k','b','r'};
cc_u={'c','m'};
units_tuned=zeros(size(allp,1),length(stims_tuned)); %initialise tuned array
units_untuned=zeros(size(allp,1),length(stims_untuned)); %initialise untuned array

for ii=1:length(stims_tuned)
    stimind=find(strcmp(stimList,stims_tuned{ii})); %your stimulus of interest, from the list in the comment above
    units_tuned(find(allp(:,stimind)<0.01),ii)=1; %all units tuned for this stimuli
end

for ii=1:length(stims_untuned)
    stimind=find(strcmp(stimList,stims_untuned{ii})); %your stimulus of interest, from the list in the comment above
    units_untuned(find(allp(:,stimind)>0.01),ii)=1; %all units untuned for this stimuli
    %units_untuned(ii).units=find(allp(:,stimind)>0.05);
end


temporal_units_original=intersect(find(sum(units_tuned,2)>1),find(sum(units_untuned,2)>0)); %find neurons that have tuning for at least 5 pitch stimuli but and at least one untuned
% temporal_units=(find(sum(units_tuned,2)>1)); %find neurons that have tuning for at least 5 pitch stimuli but and at least one untuned

temporal_units=temporal_units_original;
% new_temporal_units=zeros(length(temporal_units),1);
% new_temporal_units=[];

for t=1:length(temporal_units)
    allmrate=[];
    thisunit=temporal_units(t);
    figure(1);hold on;clf;
    bestf0=[];
    compf0=zeros(1,3);
    for kk=1:length(stims_tuned) %for every tuned stimulus
        spikerate=allDat{allpen(thisunit)}.Data(allunits(thisunit)).values; %find spike rate
        allstims=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gStim; %find stimulus
        f0=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gF0; %find frequencies
        f0=cellfun(@str2num, f0); %make F0 values numeric
        
        stim=stims_tuned{kk};
        stim2use=strcmp(allstims,stim); %find trials with this stimulus type

        trials=find(stim2use);
        f0=f0(trials);
        allf0=unique(f0); %different frequencies of all valid trials
        spikerate=spikerate(trials); %spike rate of valid trials
        mrate=[];
        serate=[];
        for jj=1:length(allf0) %for every frequency
            mrate(jj)=mean(spikerate(find(f0==allf0(jj)))); %calculate the mean spike rate
            serate(jj)=std(spikerate(find(f0==allf0(jj))))/sqrt(length(find(f0==allf0(jj)))); %calculate the SE
        end
        allmrate=[allmrate; mrate];
        [temp,bestf0ind]=max(mrate);
        bestf0(kk)=allf0(bestf0ind); %save the peak of the tuning curve here
                
        
        if isequal(stim,'alt') %& stim2use==1            
            compf0(kk-1)=(bestf0(kk-1)/2);
            if(0.74*compf0(kk-1))<bestf0(kk) && bestf0(kk)<(1.19*compf0(kk-1))
                temporal_units(t)=temporal_units(t);
            else
                temporal_units(t)=0;
            end
        end
        
%         if isequal(stim,'rand') %& stim2use==1
%             if bestf0(kk)<bestf0(kk-1)
%                 new_temporal_units(t)=temporal_units(t);
%             else
%                 disp(kk)
%             end
%         end
    end
     
[~,~,v]=find(temporal_units);
temporal_units_cleaned=v;


    
    
    
     for kk=1:length(stims_untuned)
        spikerate=allDat{allpen(thisunit)}.Data(allunits(thisunit)).values;
        allstims=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gStim;
        f0=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gF0;
        f0=cellfun(@str2num, f0); %make F0 values numeric
        
        stim=stims_untuned{kk};
        stim2use=strcmp(allstims,stim); %find trials with this stimulus type

        trials=find(stim2use);
        f0=f0(trials);
        allf0=unique(f0);
        spikerate=spikerate(trials);
        mrate=[];
        serate=[];
        for jj=1:length(allf0)
            mrate(jj)=mean(spikerate(find(f0==allf0(jj))));
            serate(jj)=std(spikerate(find(f0==allf0(jj))))/sqrt(length(find(f0==allf0(jj))));
        end
        allmrate=[allmrate; mrate];
    
     end    
end

disp(length(find(temporal_units)))
    for t=1:length(temporal_units_cleaned)
    allmrate=[];
    thisunit=temporal_units_cleaned(t);
    bestf0=[];
    for kk=1:length(stimsforaxes) %for every tuned stimulus
        spikerate=allDat{allpen(thisunit)}.Data(allunits(thisunit)).values; %find spike rate
        allstims=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gStim; %find stimulus
        f0=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gF0; %find frequencies
        f0=cellfun(@str2num, f0); %make F0 values numeric
        
        stim=stimsforaxes{kk};
        stim2use=strcmp(allstims,stim); %find trials with this stimulus type

        trials=find(stim2use);
        f0=f0(trials);
        allf0=unique(f0); %different frequencies of all valid trials
        spikerate=spikerate(trials); %spike rate of valid trials
        mrate=[];
        serate=[];
        for jj=1:length(allf0) %for every frequency
            mrate(jj)=mean(spikerate(find(f0==allf0(jj)))); %calculate the mean spike rate
            serate(jj)=std(spikerate(find(f0==allf0(jj))))/sqrt(length(find(f0==allf0(jj)))); %calculate the SE
        end
        
       
        if strcmp(stim,'low') %if the stimulus type is low harmonics
            allmrate=[allmrate; [0 0 mrate]];
        elseif strcmp(stim,'F0MaskLow')
            allmrate=[allmrate; [0 0 mrate]];
        else
            allmrate=[allmrate; mrate];
        end 
        
%         allmrate=[allmrate; mrate];
        [temp,bestf0ind]=max(mrate);
        bestf0(kk)=allf0(bestf0ind); %save the peak of the tuning curve here
    print([pwd '/Figures/' 'Tuning/TemporalNeuron_' num2str(thisunit) '.png'],'-dpng')
                
        
    end  
    
  

    figure(2);clf
    imagesc(allmrate);
    set(gca,'xtick',[1:4:length(allf0)],'xticklabels',(allf0(1:4:end)),'fontname','ariel','fontsize',40);
    stimstuneduntuned=[stims_tuned;stims_untuned];
    set(gca,'ytick',[1:length(stimsforaxes)],'yticklabels',stimsforaxes)
    colorbar      
%   title(allDat{allpen(thisunit)}.fName)
    set(gcf,'position',[679   459   592   339])
    set(gcf,'paperpositionmode','manual')
    set(gcf,'paperposition', [0 0 30 20])
    print([pwd '/Figures/' 'Tuning/TemporalNeuron_' num2str(thisunit) '_imagesc.png'],'-dpng')
    
end


    
%     pause()

%% Temporal neuron locations
fList = {'Derry_P02_Pitch' 'Derry_P03_Pitch' 'Derry_P05_Pitch' 'Derry_P08_Pitch'...
    'Ronnie_P04_Pitch' 'Ronnie_P05_Pitch' 'Ronnie_P08_Pitch' 'Ronnie_P13_Pitch'...
    'Dory_P00_Pitch' 'Dory_P01_Pitch' 'Dory_P02_Pitch' 'Dory_P04_Pitch'...
    'Noah_P01_Pitch' 'Noah_P02_Pitch' 'Noah_P03_Pitch' 'Noah_P04_Pitch' ...
    'Noah_P05_Pitch' 'Noah_P06_Pitch' 'Noah_P07_Pitch' 'Noah_P08_Pitch' };

LocPen = {'LowAAF' 'MidA1' 'LowA1' 'LowA1'...
    'LowAAF' 'LowAAF' 'LowA1' 'MidA1'...
    'LowA1' 'LowAAF' 'MidAAF' 'LowAAF'...
    'LowA1' 'LowA1' 'HighA1' 'LowA1' ...
    'LowPPF' 'LowAAF' 'LowPPF' 'High AAF'};

LocPenidx = [2 3 1 1 2 2 1 3 1 2 4 2 1 2 5 1 7 2 7 6];

%1 = LowA1, 2 = LowAAF, 3 = MidA1, 4 = MidAAF, 5 = HighA1, 6 = High AAF,
% 7 = LowPPF

TListidx = unique(LocPenidx);
TList = unique(LocPen);
finalDat = [];
clear loc_pitchUnits fname_temporalUnits
for pp=1:length(temporal_units_cleaned)
    thisunit=temporal_units_cleaned(pp);
    Dat = allDat{allpen(thisunit)};
    loc_temporalUnits{pp}=LocPen{allpen(thisunit)};
    fname_temporalUnits{pp}=fList{allpen(thisunit)};
end

fname_temporalUnits;
loc_temporalUnits; 

Tcount_all=zeros(1,length(TList));
Tcount_temporal=zeros(1,length(TList));
for c=1:length(allDat)
    numu=length(allDat{c}.ttestPval);%num units in this pen
    Tcount_all(LocPenidx(c))=Tcount_all(LocPenidx(c))+numu;
end
for a=1:length(temporal_units_cleaned)
    thisunit=temporal_units_cleaned(a);
    thisarea = LocPenidx(allpen(thisunit));
    Tcount_temporal(thisarea)=Tcount_temporal(thisarea)+1;
end
 
TList
Tcount_all
Tcount_temporal
Tcount_temporal./Tcount_all*100
%% Find harmonic encoding neurons
cd('E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData')
% load('All1wayanova.mat')
% load('allDat_allanimalsandpens.mat')
stims_tuned={'low';'CT0'}; %want tuning for these
stims_untuned={'CT40';'high'};
stimsforaxes={'low';'F0MaskLow';'high';'F0MaskHigh';'tone';'CT0';'CT5';'CT10';'CT20';'CT40';'rand';'alt'};

cc_t={'b','r'};
cc_u={'k','c'};
units_tuned=zeros(size(allp,1),length(stims_tuned)); %initialise tuned array
units_untuned=zeros(size(allp,1),length(stims_untuned)); %initialise untuned array


for ii=1:length(stims_tuned)
    stimind=find(strcmp(stimList,stims_tuned{ii})); %your stimulus of interest, from the list in the comment above
    units_tuned(find(allp(:,stimind)<0.01),ii)=1; %all units tuned for this stimuli
end
for ii=1:length(stims_untuned)
    stimind=find(strcmp(stimList,stims_untuned{ii})); %your stimulus of interest, from the list in the comment above
    units_untuned(find(allp(:,stimind)>0.01),ii)=1; %all units untuned for this stimuli
    %units_untuned(ii).units=find(allp(:,stimind)>0.05);
end

harmonic_units=intersect(find(sum(units_tuned,2)>1),find(sum(units_untuned,2)>1)); %find neurons that have tuning for at least 5 pitch stimuli but and at least one untuned

    
length(harmonic_units) %how many pitch neurons are there

for ii=1:length(harmonic_units) %for every putative neuron
    allmrate=[];
    thisunit=harmonic_units(ii);
%     figure(1);hold on;clf
    bestf0=[];
    for kk=1:length(stimsforaxes) %for every tuned stimulus
        spikerate=allDat{allpen(thisunit)}.Data(allunits(thisunit)).values; %find spike rate
        allstims=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gStim; %find stimulus
        f0=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gF0; %find frequencies
        f0=cellfun(@str2num, f0); %make F0 values numeric
        
        stim=stimsforaxes{kk};
        stim2use=strcmp(allstims,stim); %find trials with this stimulus type

        trials=find(stim2use);
        f0=f0(trials);
        allf0=unique(f0); %different frequencies of all valid trials
        spikerate=spikerate(trials); %spike rate of valid trials
        mrate=[];
        serate=[];
        for jj=1:length(allf0) %for every frequency
            mrate(jj)=mean(spikerate(find(f0==allf0(jj)))); %calculate the mean spike rate
            serate(jj)=std(spikerate(find(f0==allf0(jj))))/sqrt(length(find(f0==allf0(jj)))); %calculate the SE
        end
        
       
        if strcmp(stim,'low') %if the stimulus type is low harmonics
            allmrate=[allmrate; [0 0 mrate]];
        elseif strcmp(stim,'F0MaskLow')
            allmrate=[allmrate; [0 0 mrate]];
        else
            allmrate=[allmrate; mrate];
        end 
        
%         allmrate=[allmrate; mrate];
        [temp,bestf0ind]=max(mrate);
        bestf0(kk)=allf0(bestf0ind); %save the peak of the tuning curve here
                
        
    end  
    
  

    figure(2);clf
    imagesc(allmrate);
    set(gca,'xtick',[1:4:length(allf0)],'xticklabels',(allf0(1:4:end)),'fontname','ariel','fontsize',40);
    stimstuneduntuned=[stims_tuned;stims_untuned];
    set(gca,'ytick',[1:length(stimsforaxes)],'yticklabels',stimsforaxes)
    colorbar      
%   title(allDat{allpen(thisunit)}.fName)
    set(gcf,'position',[679   459   592   339])
    set(gcf,'paperpositionmode','manual')
    set(gcf,'paperposition', [0 0 30 20])
%     print([pwd '/Figures/' 'Tuning/HarmonicNeuron_' num2str(thisunit) '_imagesc.png'],'-dpng')
    


    
end
%% Harmonic neuron locations
fList = {'Derry_P02_Pitch' 'Derry_P03_Pitch' 'Derry_P05_Pitch' 'Derry_P08_Pitch'...
    'Ronnie_P04_Pitch' 'Ronnie_P05_Pitch' 'Ronnie_P08_Pitch' 'Ronnie_P13_Pitch'...
    'Dory_P00_Pitch' 'Dory_P01_Pitch' 'Dory_P02_Pitch' 'Dory_P04_Pitch'...
    'Noah_P01_Pitch' 'Noah_P02_Pitch' 'Noah_P03_Pitch' 'Noah_P04_Pitch' ...
    'Noah_P05_Pitch' 'Noah_P06_Pitch' 'Noah_P07_Pitch' 'Noah_P08_Pitch' };

LocPen = {'LowAAF' 'MidA1' 'LowA1' 'LowA1'...
    'LowAAF' 'LowAAF' 'LowA1' 'MidA1'...
    'LowA1' 'LowAAF' 'MidAAF' 'LowAAF'...
    'LowA1' 'LowA1' 'HighA1' 'LowA1' ...
    'LowPPF' 'LowAAF' 'LowPPF' 'High AAF'};

LocPenidx = [2 3 1 1 2 2 1 3 1 2 4 2 1 2 5 1 7 2 7 6];

%1 = LowA1, 2 = LowAAF, 3 = MidA1, 4 = MidAAF, 5 = HighA1, 6 = High AAF,
% 7 = LowPPF

HListidx = unique(LocPenidx);
HList = unique(LocPen);
finalDat = [];
clear loc_harmonicUnits fname_harmonicUnits
for pp=1:length(harmonic_units)
    thisunit=harmonic_units(pp);
    Dat = allDat{allpen(thisunit)};
    loc_harmonicUnits{pp}=LocPen{allpen(thisunit)};
    fname_harmonicUnits{pp}=fList{allpen(thisunit)};
end

fname_harmonicUnits
loc_harmonicUnits 

Hcount_all=zeros(1,length(HList));
Pcount_harmonic=zeros(1,length(HList));
for pp=1:length(allDat)
    numu=length(allDat{pp}.ttestPval);%num units in this pen
    Hcount_all(LocPenidx(pp))=Hcount_all(LocPenidx(pp))+numu;
end
for pp=1:length(harmonic_units)
    thisunit=harmonic_units(pp);
    thisarea = LocPenidx(allpen(thisunit));
    Hcount_harmonic(thisarea)=Hcount_harmonic(thisarea)+1;
end
 
HList
Hcount_all
Hcount_pitch
Hcount_pitch./Hcount_all*100
%% Harmonic and temporal neuron percentages
%Run blocks to find each first
cd ('E:\Dan analysis\Daniel_ephys_analysis\pitch_ephys\DansMATLABData')
% load('All1wayanova.mat');
stims={'high'; 'alt'}; %your stimulus of interest, from the list in the comment above
percenttuned=zeros(1,2);
% for ii=1:length(stims) %for every stimulus
%     stimind=find(strcmp(stimList,stims{ii})); %find the stimulus in the list
%     fsig=find(allp(:,stimind)<0.05); %find all the units tuned to this stimulus
%     percenttuned(ii)=length(fsig)/size(allp,1)*100; %find the percentage tuned from these
% end
percenttuned(1)=length(temporal_units_cleaned)/size(allp,1)*100;
percenttuned(2)=length(harmonic_units)/size(allp,1)*100;

figure(2);clf;hold on; %make a figure window
bar([1:length(stims)]*4,percenttuned) %make a bar chart for all stim on x, against % tuned
bar([1:2]*4,percenttuned) %make a bar chart for all stim on x, against % tuned

set(gca,'xtick',[1:length(stims)]*4,'xticklabels',stims,'fontname','ariel','fontsize',40);
%xlabel('F0 (Hz)','fontname','ariel','fontsize',32);
ylabel('% of F0-sensitive neurons','fontname','ariel','fontsize',40);
box off

% set(gcf,'position',[582   503   781   279])
% set(gcf,'paperpositionmode','manual')
% set(gcf,'paperposition',[0 0 30 25])
% print([pwd '/Figures/' 'Tuning/PercentTuned_highalt.png'],'-dpng')

stims2={''; 'alt'}; 
controlstim='high'; 
stimind2=find(strcmp(stimList,controlstim)); 
fsig_control=find(allp(:,stimind2)<0.05); %indices of neurons tuned to F0 for the control stim
percenttuned_both=[];
% for ii=1:length(stims2)
%     stimind3=find(strcmp(stimList,stims2{ii})); %your stimulus of interest, from the list in the comment above
%     fsig3=find(allp(:,stimind3)<0.05); %index of all neurons tuned to F0 for other stims
%     fsigboth=intersect(fsig3,fsig_control); %all units tuned to both
%     percenttuned_both(ii)=length(fsigboth)/size(allp,1)*100; 
% end
% bar([1:length(stims2)]*4,percenttuned_both,'r')
% xlim([2 length(stims2)*4+6])
hold off
%% Alt tuning percentages

stim='high';
controlstim='alt';


stimind=find(strcmp(stimList,stim)); %find the stimulus in the list
fsig=find(allp(:,stimind)<0.05); %find all the units tuned to this stimulus
notfsig=find(allp(:,stimind)>=0.05); %find all the units tuned to this stimulus


similaralt=0;
loweralt=0;
higheralt=0;

for n=1:length(fsig)
    allmrate=[];
    thisunit=fsig(n);
%     figure(1);hold on;clf;
    bestf0=[];
    spikerate=allDat{allpen(thisunit)}.Data(allunits(thisunit)).values; %find spike rate
    allstims=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gStim; %find stimulus
    f0=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gF0; %find frequencies
    f0=cellfun(@str2num, f0); %make F0 values numeric

    stim2use=strcmp(allstims,stim); %find trials with this stimulus type

    trials=find(stim2use);
    f0=f0(trials);
    allf0=unique(f0); %different frequencies of all valid trials
    spikerate=spikerate(trials); %spike rate of valid trials
    mrate=[];
    serate=[];
    for jj=1:length(allf0) %for every frequency
        mrate(jj)=mean(spikerate(find(f0==allf0(jj)))); %calculate the mean spike rate
        serate(jj)=std(spikerate(find(f0==allf0(jj))))/sqrt(length(find(f0==allf0(jj)))); %calculate the SE
    end
    allmrate=[allmrate; mrate];
    [temp,bestf0ind]=max(mrate);
    bestf0(n)=allf0(bestf0ind); %save the peak of the tuning curve here

        

    
 
    spikerate=allDat{allpen(thisunit)}.Data(allunits(thisunit)).values; %find spike rate
    allstims=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gStim; %find stimulus
    f0=allDat{allpen(thisunit)}.Data(allunits(thisunit)).gF0; %find frequencies
    f0=cellfun(@str2num, f0); %make F0 values numeric

    stim2use=strcmp(allstims,controlstim); %find trials with this stimulus type

    trials=find(stim2use);
    f0=f0(trials);
    allf0=unique(f0); %different frequencies of all valid trials
    spikerate=spikerate(trials); %spike rate of valid trials
    mrate=[];
    serate=[];
    for jj=1:length(allf0) %for every frequency
        mrate(jj)=mean(spikerate(find(f0==allf0(jj)))); %calculate the mean spike rate
        serate(jj)=std(spikerate(find(f0==allf0(jj))))/sqrt(length(find(f0==allf0(jj)))); %calculate the SE
    end
    allmrate=[allmrate; mrate];
    [temp,bestf0ind]=max(mrate);
    controlbestf0(n)=allf0(bestf0ind); %save the peak of the tuning curve here
    
    
    
    
    if 0.84*controlbestf0(n)<bestf0(n) & bestf0(n)<1.19*controlbestf0(n)
        similaralt=similaralt+1;
    elseif controlbestf0(n)<bestf0(n)
        loweralt=loweralt+1;
    elseif controlbestf0(n)>bestf0(n)
        higheralt=higheralt+1;
    end
    

    constimind=find(strcmp(stimList,controlstim)); %find the stimulus in the list
    fsigcontrol=find(allp(:,constimind)<0.05);
    notfsigcontrol=find(allp(:,constimind)>=0.05); %find all the units tuned to this stimulus
    
    
    
    highnotalt=intersect(fsig,notfsigcontrol);
    altnothigh=intersect(fsigcontrol,notfsig);
    
end

eitheralthigh=unique(cat(1,fsigcontrol,fsig));

highTotal=length(fsig) %high responsive total
similaralt
loweralt %these 3 should add up to above
higheralt
length(fsigcontrol) %alt responsive total
length(highnotalt) %just high
length(altnothigh) %just alt
length(intersect(fsig,fsigcontrol)) %both
length(eitheralthigh)


percentclose=(similaralt/highTotal)*100;
percentlower=(loweralt/highTotal)*100;
percenthigher=(higheralt/highTotal)*100;
percentnottuned=(length(highnotalt)/highTotal)*100;

toPlot=[percentclose;percentlower;percenthigher;percentnottuned];

bar(toPlot)



%% Before this, run unit selection, then pitch neurons, then location of pitch neurons
%Works, is the same as below I just fixed them both without realising they
%were the same - Dan
unitMat = {};
c=1;
stimNames1 = {'low', 'F0MaskLow'};
stimNames2 = {'high', 'F0MaskHigh', 'alt', 'rand', 'allHarm'};
h1 = figure; h2 = figure; h3 = figure; h4 = figure;
for imp = 1:length(allDat) %for all data
    Dat = allDat{imp}; %get data
    tmp = regexp(Dat.fName,'_','split'); %split the fName into the 3 constituents which were separated by '_'
    an = tmp{1}; %get animal name
    pen = tmp{2}; %get penetration number (form = Pn where n = 2 digit number)
    pLoc = LocPen(strcmp(fList,Dat.fName)); %find the location of the penetration
    nU = length(Dat.idxResp); 
    
    for i = 1:nU, %for every responsive neuron
        
        for j = 1:length(stimNames1), 
            [meanMatLow(j,:), errMatLow(j,:), FListLow, pvalLow(j), CFlow(j)] = calc_TC(Dat.Data(i),gStim,gF0,stimNames1{j}); %make arrays for these variables
        end
        meanMatLow2 = [nan(2,2) meanMatLow]; %adds 2x2 NaN to start of array
        errMatLow2 = [nan(2,2) errMatLow];
        for j = 1:length(stimNames2),
            [meanMat(j,:) errMat(j,:), FList, pval(j), CF(j)] = calc_TC(Dat.Data(i),gStim,gF0,stimNames2{j}); %same as above for the other stimuli
        end
%         meanMat2 = [meanMat([6 8 9 10 7 2],:);meanMatLow2(2,:); meanMat([5 3 1],:); meanMatLow2(1,:); meanMat([4 11],:)]; 
%         errMat2 = [errMat([6 8 9 10 7 2],:); errMatLow2(2,:); errMat([5 3 1],:); errMatLow2(1,:); errMat([4 11],:)];
%         BF0 = [CF([6 8 9 10 7 2]) CFlow(2) CF([5 3 1]) CFlow(1) CF([4 11])];

%the above is broken, but appears to have no function so I've commented it
%out along with unitMat below - Dan

        an1pVal = Dat.DataAnova1(i).p;
        resp = Dat.idxRespTtest(i);
%         unitMat = [unitMat; [{an} {pen} {pLoc} {i} {resp} {an1pVal} {meanMat2} {errMat2} {BF0}]];
        c = c+1;
        
           figure(h1);
    errorbar(FListLow,meanMatLow(1,:), errMatLow(1,:))
    hold on;
    errorbar(FListLow,meanMatLow(2,:), errMatLow(2,:))
    hold off;
    legend(stimNames1);
    
    figure(h2);
    errorbar(FListHigh,meanMatHigh(1,:), errMatHigh(1,:))
    hold on;
    errorbar(FListHigh,meanMatHigh(2,:), errMatHigh(2,:))
    hold off;
    legend(stimNames2([1 2]));
    
    % Phase modification effect
    figure(h3)
    errorbar(FListHigh,meanMatHigh(1,:), errMatHigh(1,:))
    hold on;
    errorbar(FListHigh,meanMatHigh(3,:), errMatHigh(3,:))
    errorbar(FListHigh,meanMatHigh(4,:), errMatHigh(4,:))
    hold off;
    legend(stimNames2([1 3 4]));
    
    % All harm vs. high harm vs. low harm
    figure(h4)
    errorbar(FListHigh,meanMatHigh(5,:), errMatHigh(5,:))
    hold on;
    errorbar(FListHigh,meanMatHigh(1,:), errMatHigh(1,:))
    errorbar(FListLow,meanMatLow(1,:), errMatLow(1,:))
    hold off;
    legend([stimNames2([5 1]) stimNames1(1)])
    
    pause();
    end
   
end

 %% Plot Tuning curve, raster and CF for stimuli with and without masker
% Run Unit Selection first
%works

for r = 1:length(allDat)
    uIdx = find(allDat{r}.idxResp);
end


%uIdx = find(idxResp);
CFValues = zeros(length(uIdx),2);

stimNames1 = {'low', 'F0MaskLow'};
stimNames2 = {'high', 'F0MaskHigh', 'alt', 'rand', 'allHarm'};

h1 = figure; h2 = figure; h3 = figure; h4 = figure;

for imp = 1:length(allDat) %for all data
    Dat = allDat{imp}; %get data
    tmp = regexp(Dat.fName,'_','split'); %split the fName into the 3 constituents which were separated by '_'
    an = tmp{1}; %get animal name
    pen = tmp{2}; %get penetration number (form = Pn where n = 2 digit number)
    pLoc = LocPen(strcmp(fList,Dat.fName)); %find the location of the penetration
    nU = length(Dat.idxResp); 

for i = 1:length(nU),
    
    for j = 1:length(stimNames1),
        [meanMatLow(j,:), errMatLow(j,:), FListLow] = calc_TC(Dat.Data(uIdx(i)),gStim,gF0,stimNames1{j});
    end
    
    for j = 1:length(stimNames2),
        [meanMatHigh(j,:), errMatHigh(j,:), FListHigh] = calc_TC(Dat.Data(uIdx(i)),gStim,gF0,stimNames2{j});
    end
    
    % Masking effect
    figure(h1);
    errorbar(FListLow,meanMatLow(1,:), errMatLow(1,:))
    hold on;
    errorbar(FListLow,meanMatLow(2,:), errMatLow(2,:))
    hold off;
    legend(stimNames1);
    
    figure(h2);
    errorbar(FListHigh,meanMatHigh(1,:), errMatHigh(1,:))
    hold on;
    errorbar(FListHigh,meanMatHigh(2,:), errMatHigh(2,:))
    hold off;
    legend(stimNames2([1 2]));
    
    % Phase modification effect
    figure(h3)
    errorbar(FListHigh,meanMatHigh(1,:), errMatHigh(1,:))
    hold on;
    errorbar(FListHigh,meanMatHigh(3,:), errMatHigh(3,:))
    errorbar(FListHigh,meanMatHigh(4,:), errMatHigh(4,:))
    hold off;
    legend(stimNames2([1 3 4]));
    
    % All harm vs. high harm vs. low harm
    figure(h4)
    errorbar(FListHigh,meanMatHigh(5,:), errMatHigh(5,:))
    hold on;
    errorbar(FListHigh,meanMatHigh(1,:), errMatHigh(1,:))
    errorbar(FListLow,meanMatLow(1,:), errMatLow(1,:))
    hold off;
    legend([stimNames2([5 1]) stimNames1(1)])
    
    pause();
    end
end
%% Plot mimicking the behavior results
% Run Unit Selection first
%works

%uIdx = find(idxResp);


for r = 1:length(allDat)
    uIdx = find(allDat{r}.idxResp);
end

CFValues = zeros(length(uIdx),2);

stimNames1 = {'low', 'F0MaskLow'};
stimNames2 = {'high', 'F0MaskHigh', 'alt', 'rand', 'allHarm'};

a = 0.05;

finalMat = [];
CFMat = [];
finalMatNames = {'Standard' 'Resolved' 'Resolved Mskd' 'Unresolved' 'Unresolved Mskd' 'Unresolved alt' 'Unresolved rand'};
nU = length(uIdx);
for i = 1:nU,
    
    for j = 1:length(stimNames1),
        [meanMatLow(j,:), errMatLow(j,:), FListLow, pvalLow(j,:), CFlow(j)] = calc_TC(Dat.Data(uIdx(i)),gStim,gF0,stimNames1{j});
    end
    
    for j = 1:length(stimNames2),
        [meanMat(j,:), errMat(j,:), FList, pval(j,:), CF(j)] = calc_TC(Dat.Data(uIdx(i)),gStim,gF0,stimNames2{j});
    end
    
    [c,m,h,gnames] = multcompare(Dat.Data(uIdx(i)).stats,'Alpha',0.05,'CType','tukey-kramer','Display','off','Dimension',[2]);
    
    finalMat = [finalMat ; [(sum(pval(5,:)<a)>0) (sum(pvalLow(1,:)<a)>0) (sum(pvalLow(2,:)<a)>0) ...
        (sum(pval(1,:)<a)>0) (sum(pval(2,:)<a)>0) (sum(pval(3,:)<a)>0) (sum(pval(4,:)<a)>0)]];
    CFMat = [CFMat; [CF(5) CFlow(1:2) CF(1:4)]];
end

finalMat2 = finalMat(logical(finalMat(:,1)),:);
CFMat = CFMat(logical(finalMat(:,1)),:);

figure;
bar(sum(finalMat2)./size(finalMat2,1)*100)
set(gca,'xticklabel',finalMatNames,'xticklabelRotation',45)
ylabel('% of units')

CFMat = CFMat(logical(finalMat2(:,4)) & logical(finalMat2(:,6)) & logical(finalMat2(:,7)),:);
figure;
boxplot(CFMat(:,[1 4 6 7]))
ylabel('CF (Hz)');
set(gca,'xticklabel',finalMatNames([1 4 6 7]),'xticklabelRotation',45)

% Plot difference in prefered F0 in octave
idx = [4 6 7];
for i = 1:3,
%     d = log2(CFMat(:,1)./CFMat(:,idx(i)));
    figure;
%     hist(d,10)
    scatterQ(CFMat(:,1),CFMat(:,idx(i)),1)
end

%% Plot mimicking the behavior results - Penetration by penetration
% Run Unit Selection first
%works, but what is it plotting and why is Noah not there (or Dory)

stimNames1 = {'low', 'F0MaskLow'};
stimNames2 = {'high', 'F0MaskHigh', 'alt', 'rand', 'allHarm'};

a = 0.05;

finalMat = {};
CFMat = {};
finalMatNames = {'Standard' 'Resolved' 'Resolved Mskd' 'Unresolved' 'Unresolved Mskd' 'Unresolved alt' 'Unresolved rand'};

pairList = [8 13; 8 1; 8 11; 8 7; 8 10 ; 8 6;8 9; 8 12]; %what is this

for imp = 1:length(allDat),
    Dat = allDat{imp};
    uIdx = find(Dat.idxResp);
    CFValues = zeros(length(uIdx),2);
    nU = length(uIdx);
    nAllU = length(Dat.idxResp);
    
    finalMat = zeros(nU,7);
    finalMatAH = zeros(nU,7);
    CFMat = zeros(nU,7);
    
    for i = 1:nU,
        
        for j = 1:length(stimNames1),
            [meanMatLow(j,:), errMatLow(j,:), FListLow, pvalLow(j,:), CFlow(j)] = calc_TC(Dat.Data(uIdx(i)),gStim,gF0,stimNames1{j});
        end
        
        for j = 1:length(stimNames2),
            [meanMat(j,:), errMat(j,:), FList, pval(j,:), CF(j)] = calc_TC(Dat.Data(uIdx(i)),gStim,gF0,stimNames2{j});
        end
        
%         [c,m,h,gnames] = multcompare(Dat.Data(uIdx(i)).stats,'Alpha',0.05,'CType','tukey-kramer','Display','on','Dimension',[1]);
        [c,m,h,gnames] = multcompare(Dat.Data(uIdx(i)).anova2.stats,'Alpha',0.05,'CType','tukey-kramer','Display','off','Estimate','row');
        hC = c(:,3) > 0 | c(:,5)<0;
        for pl = 1:size(pairList,1),
            finalMatAH(i,pl) = hC(((c(:,1) == pairList(pl,1)) | (c(:,2) == pairList(pl,1))) & ((c(:,1) == pairList(pl,2)) | (c(:,2) == pairList(pl,2))));
        end
                
        finalMat(i,:) = [(sum(pval(5,:)<a)>0) (sum(pvalLow(1,:)<a)>0) (sum(pvalLow(2,:)<a)>0) ...
            (sum(pval(1,:)<a)>0) (sum(pval(2,:)<a)>0) (sum(pval(3,:)<a)>0) (sum(pval(4,:)<a)>0)];
        CFMat(i,:) = [CF(5) CFlow(1:2) CF(1:4)];
    end
    allMat{imp}.finalMat = finalMat;
    allMat{imp}.finalMatAH = finalMatAH;
    allMat{imp}.CFMat = CFMat;
    allMat{imp}.nU = nU;
    allMat{imp}.nAllU = nAllU;
end

allCFMat = [];
cellProp = [];
for imp = 1:length(allDat),
    
    idx = logical(allMat{imp}.finalMat(:,1));
    finalMat2 = allMat{imp}.finalMat(logical(allMat{imp}.finalMat(:,1)),:);
    
%     finalMat2 = allMat{imp}.finalMatAH(idx,:);
%     finalMat2 = double(~logical(finalMat2));

    CFMat = allMat{imp}.CFMat(idx,:);
    
    cellProp(imp,:) = sum(finalMat2)./size(finalMat2,1)*100;
    
    % Plot each penetration
    %     figure;
    %     bar(sum(finalMat2)./size(finalMat2,1)*100)
    %     set(gca,'xticklabel',finalMatNames,'xticklabelRotation',45)
    %     ylabel('% of units')
    
    
    CFMat = CFMat(logical(finalMat2(:,4)) & logical(finalMat2(:,6)) & logical(finalMat2(:,7)),:);
    
    % Plot each penetration
    % igure;
    % boxplot(CFMat(:,[1 4 6 7]))
    % ylabel('CF (Hz)');
    % set(gca,'xticklabel',finalMatNames([1 4 6 7]),'xticklabelRotation',45)
    % 
    % % Plot difference in prefered F0 in octave
    % idx = [4 6 7];
    % for i = 1:3,
    % %     d = log2(CFMat(:,1)./CFMat(:,idx(i)));
    %     figure;
    % %     hist(d,10)
    %     scatterQ(CFMat(:,1),CFMat(:,idx(i)),1);
    % end

    allCFMat = [allCFMat; CFMat];
end
cellProp = cellProp(~isnan(cellProp(:,1)),:);
% Plot average over all implantations
figure;
errorbar(1:7,mean(cellProp),ste(cellProp)); %changed errorbarq to errorbar
set(gca,'xticklabel',finalMatNames,'xticklabelRotation',45)
ylabel('% of units')


figure;
boxplot(allCFMat(:,[1 4 6 7]))
ylabel('CF (Hz)');
set(gca,'xticklabel',finalMatNames([1 4 6 7]),'xticklabelRotation',45)

% Plot difference in prefered F0 in octave
idx = [4 6 7];
for i = 1:3,
%     d = log2(CFMat(:,1)./CFMat(:,idx(i)));
    figure;
%     hist(d,10)
    scatterQ(allCFMat(:,1),allCFMat(:,idx(i)),1);
end

for i = 1:6,
    fprintf('%s: %2.2f%% Selective out of %d units\n',fList{i},allMat{i}.nU/allMat{i}.nAllU*100,allMat{i}.nAllU);
end

%% Plot mimicking the behavior results - Pooled all units
% Run Unit Selection first
%same questions as above
stimNames1 = {'low', 'F0MaskLow'};
stimNames2 = {'high', 'F0MaskHigh', 'alt', 'rand', 'allHarm'};

a = 0.05;

finalMat = {};
CFMat = {};
finalMatNames = {'All Harm.' 'Resolved' 'Resolved Mskd' 'Unresolved' 'Unresolved Mskd' 'Unresolved alt' 'Unresolved rand'};

pairList = [8 13; 8 1; 8 11; 8 7; 8 10 ; 8 6;8 9; 8 12];

for imp = 1:length(allDat),
    Dat = allDat{imp};
    uIdx = find(Dat.idxResp);
    CFValues = zeros(length(uIdx),2);
    nU = length(uIdx);
    nAllU = length(Dat.idxResp);
    
    finalMat = zeros(nU,7);
    finalMatAH = zeros(nU,7);
    CFMat = zeros(nU,7);
    
    for i = 1:nU,
        
        for j = 1:length(stimNames1),
            [meanMatLow(j,:), errMatLow(j,:), FListLow, pvalLow(j,:), CFlow(j)] = calc_TC(Dat.Data(uIdx(i)),gStim,gF0,stimNames1{j});
        end
        
        for j = 1:length(stimNames2),
            [meanMat(j,:), errMat(j,:), FList, pval(j,:), CF(j)] = calc_TC(Dat.Data(uIdx(i)),gStim,gF0,stimNames2{j});
        end
        
%         [c,m,h,gnames] = multcompare(Dat.Data(uIdx(i)).stats,'Alpha',0.05,'CType','tukey-kramer','Display','on','Dimension',[1]);
        [c,m,h,gnames] = multcompare(Dat.Data(uIdx(i)).anova2.stats,'Alpha',0.05,'CType','tukey-kramer','Display','off','Estimate','row');
        hC = c(:,3) > 0 | c(:,5)<0;
        for pl = 1:size(pairList,1),
            finalMatAH(i,pl) = hC(((c(:,1) == pairList(pl,1)) | (c(:,2) == pairList(pl,1))) & ((c(:,1) == pairList(pl,2)) | (c(:,2) == pairList(pl,2))));
        end
                
        finalMat(i,:) = [(sum(pval(5,:)<a)>0) (sum(pvalLow(1,:)<a)>0) (sum(pvalLow(2,:)<a)>0) ...
            (sum(pval(1,:)<a)>0) (sum(pval(2,:)<a)>0) (sum(pval(3,:)<a)>0) (sum(pval(4,:)<a)>0)];
        CFMat(i,:) = [CF(5) CFlow(1:2) CF(1:4)];
    end
    allMat{imp}.finalMat = finalMat;
    allMat{imp}.finalMatAH = finalMatAH;
    allMat{imp}.CFMat = CFMat;
    allMat{imp}.nU = nU;
    allMat{imp}.nAllU = nAllU;
end

allCFMat = [];
cellProp = [];
finalMat2 = [];
for imp = 1:length(allDat),
    
%     idx = logical(allMat{imp}.finalMat(:,1));
    idx = true(size(allMat{imp}.finalMat(:,1),1),1);
    finalMat2 = [finalMat2; allMat{imp}.finalMat(idx,:)];

    CFMat = allMat{imp}.CFMat(idx,:);
    allCFMat = [allCFMat; CFMat];
    
%     cellProp(imp,:) = sum(finalMat2)./size(finalMat2,1)*100
    
    % Plot each penetration
    %     figure;
    %     bar(sum(finalMat2)./size(finalMat2,1)*100)
    %     set(gca,'xticklabel',finalMatNames,'xticklabelRotation',45)
    %     ylabel('% of units')
    
    
%     CFMat = CFMat(logical(finalMat2(:,4)) & logical(finalMat2(:,6)) & logical(finalMat2(:,7)),:);
    
    % Plot each penetration
    % igure;
    % boxplot(CFMat(:,[1 4 6 7]))
    % ylabel('CF (Hz)');
    % set(gca,'xticklabel',finalMatNames([1 4 6 7]),'xticklabelRotation',45)
    % 
    % % Plot difference in prefered F0 in octave
    % idx = [4 6 7];
    % for i = 1:3,
    % %     d = log2(CFMat(:,1)./CFMat(:,idx(i)));
    %     figure;
    % %     hist(d,10)
    %     scatterQ(CFMat(:,1),CFMat(:,idx(i)),1);
    % end

    
end
cellProp = sum(finalMat2)./size(finalMat2,1)*100;
% Plot average over all implantations
figure;
bar(1:4,cellProp([3 5:7]));
set(gca,'xticklabel',finalMatNames([3 5:7]),'xticklabelRotation',45)
ylabel('% of pitch units')
set(gca,'FontWeight','bold','YGrid','on')


figure('position',[680   364   560   614]);
boxplot(allCFMat(:,[1 4 6 7]))
ylabel('Pref. F0 (Hz)');
set(gca,'xticklabel',finalMatNames([1 4 6 7]),'xticklabelRotation',45)
set(gca,'FontSize',20,'FontWeight','bold','ygrid','on')

% Plot difference in prefered F0 in octave
idx = [4 6 7];
for i = 1:3,
%     d = log2(CFMat(:,1)./CFMat(:,idx(i)));
    figure;
%     hist(d,10)
    scatterQ(allCFMat(:,1),allCFMat(:,idx(i)),1);
    set(gca,'FontSize',20,'FontWeight','bold')
    xlabel(['Pref. F0 - ' finalMatNames{1}]);
    ylabel(['Pref. F0 - ' finalMatNames{idx(i)}]);
    [H,p] = ttest(allCFMat(:,1),allCFMat(:,idx(i)));
    fprintf('p = %d\n',p)
end

for i = 1:6,
    fprintf('%s: %2.2f%% Selective out of %d units\n',fList{i},allMat{i}.nU/allMat{i}.nAllU*100,allMat{i}.nAllU);
end

%% All Harm vs. tones - Pooled all units
% Run Unit Selection first
%not working
stimNames1 = {'allHarm', 'tone'};
% stimNames2 = {'high', 'F0MaskHigh', 'alt', 'rand', 'allHarm'};

a = 0.05;

finalMat = {};
CFMat = {};
finalMatNames = {'All Harm.' 'Tones'};

% pairList = [8 13; 8 1; 8 11; 8 7; 8 10 ; 8 6;8 9; 8 12];
h = figure;
for imp = 1:length(allDat),
    Dat = allDat{imp};
    uIdx = find(Dat.idxResp);
    CFValues = zeros(length(uIdx),2); 
    nU = length(uIdx);
    nAllU = length(Dat.idxResp);
    
    finalMat = zeros(nU,2); 
%     finalMatAH = zeros(nU,7);
    CFMat = zeros(nU,2); 
    
    for i = 1:nU,
        
        for j = 1:length(stimNames1),
            [meanMat(j,:), errMat(j,:), FList, pval(j,:), CF(j)] = calc_TC(Dat.Data(uIdx(i)),gStim,gF0,stimNames1{j});
        end
        
%         [c,m,h,gnames] = multcompare(Dat.Data(uIdx(i)).anova2.stats,'Alpha',0.05,'CType','tukey-kramer','Display','off','Estimate','row');
%         hC = c(:,3) > 0 | c(:,5)<0;
%         for pl = 1:size(pairList,1),
%             finalMatAH(i,pl) = hC(((c(:,1) == pairList(pl,1)) | (c(:,2) == pairList(pl,1))) & ((c(:,1) == pairList(pl,2)) | (c(:,2) == pairList(pl,2))));
%         end
                
        finalMat(i,:) = [(sum(pval(1,:)<a)>0) (sum(pval(2,:)<a)>0)];
        CFMat(i,:) = CF; 
        
        if (sum(pval(1,:)<a)>0) && (sum(pval(2,:)<a)>0)
            figure(h)
            errorbar(FList,meanMat(1,:),errMat(1,:));
            hold on
            errorbar(FList,meanMat(2,:),errMat(2,:));
            hold off
            set(gca,'XScale','log')
            pause();
        end
    end
    allMat{imp}.finalMat = finalMat;
%     allMat{imp}.finalMatAH = finalMatAH;
    allMat{imp}.CFMat = CFMat;
    allMat{imp}.nU = nU;
    allMat{imp}.nAllU = nAllU;
end

allCFMat = [];
cellProp = [];
finalMat2 = [];
for imp = 1:length(allDat),
    
    idx = logical(allMat{imp}.finalMat(:,1)) & logical(allMat{imp}.finalMat(:,2));
%     idx = true(size(allMat{imp}.finalMat(:,1),1),1);
    finalMat2 = [finalMat2; allMat{imp}.finalMat(idx,:)];

    CFMat = allMat{imp}.CFMat(idx,:);
    allCFMat = [allCFMat; CFMat];

end
cellProp = sum(finalMat2)./size(finalMat2,1)*100;
% Plot average over all implantations
figure;
bar(1:2,cellProp);
set(gca,'xticklabel',finalMatNames,'xticklabelRotation',45)
ylabel('% of units')


figure;
boxplot(allCFMat(:,[1 2]))
ylabel('Pref. F0 (Hz)');
set(gca,'xticklabel',finalMatNames([1 2]),'xticklabelRotation',45)
set(gca,'FontSize',20,'FontWeight','bold')

% Plot difference in preferred F0 in octave
    figure;
%     hist(d,10)
    scatterQ(allCFMat(:,1),allCFMat(:,2),1);
    set(gca,'FontSize',20,'FontWeight','bold')
    xlabel(['Pref. F0 - ' finalMatNames{1}]);
    ylabel(['Pref. F0 - ' finalMatNames{2}]);
    [H,p] = ttest(allCFMat(:,1),allCFMat(:,2));
    fprintf('p = %d\n',p)

%% Click train data - all units pooled
%works

stimNames1 = {'allHarm', 'CT0', 'CT5', 'CT10', 'CT20', 'CT40'};

a = 0.05;

col = lines;

finalMat = {};
CFMat = {};
finalMatNames = {'All Harm.', 'no jit.', '5% jit.', '10% jit.', '20% jit.', '40% jit.'};
h = figure;
for imp = 1:length(allDat),
    Dat = allDat{imp};
    uIdx = find(Dat.idxResp);
    CFValues = zeros(length(uIdx),2);
    nU = length(uIdx);
    nAllU = length(Dat.idxResp);
    
    finalMat = zeros(nU,6);
    CFMat = zeros(nU,6);
    
    for i = 1:nU,
        
        for j = 1:length(stimNames1),
            [meanMat(j,:), errMat(j,:), FList, pval(j,:), CF(j)] = calc_TC(Dat.Data(uIdx(i)),gStim,gF0,stimNames1{j});
        end
               
        finalMat(i,:) = pval < a;
        CFMat(i,:) = CF;
        
%         if logical(finalMat(i,2)) && logical(finalMat(i,3)),
%         figure(h); % Browsing for examples
%         clf
%         hold on;
%         for kk = 2:6,
% %             errorbar(FList,meanMat(kk,:),errMat(kk,:));
% % error_area( x,y,e,curve_col,area_col,linewidth,scale,X_ticks_label, Y_ticks_label)
%             error_area(FList,meanMat(kk,:),errMat(kk,:),col(kk-1,:),col(kk-1,:),1,'xlog',[200 400 800 1600 3200 6400]);
%         end
%         hold off;
% %         legend(finalMatNames(2:6));
% %         set(gca,'XScale','log')
% %         xlim([200 5000])
% %         set(gca,'XTick',[200 400 800 1600 3200])
%         xlabel('F0 (Hz)')
%         ylabel('Firing rate (spike/sec)')
%         set(gca,'FontSize',20,'FontWeight','bold')
%         pause();
%         end
    end
    allMat{imp}.finalMat = finalMat;
    allMat{imp}.CFMat = CFMat;
    allMat{imp}.nU = nU;
    allMat{imp}.nAllU = nAllU;
end

allCFMat = [];
cellProp = [];
finalMat2 = [];
for imp = 1:length(allDat),
    
    idx = logical(allMat{imp}.finalMat(:,1));
%     idx = true(size(allMat{imp}.finalMat(:,1),1),1);
    finalMat2 = [finalMat2; allMat{imp}.finalMat(idx,:)];

    CFMat = allMat{imp}.CFMat(idx,:);
    allCFMat = [allCFMat; CFMat];
    
%     cellProp(imp,:) = sum(finalMat2)./size(finalMat2,1)*100
    
    % Plot each penetration
    %     figure;
    %     bar(sum(finalMat2)./size(finalMat2,1)*100)
    %     set(gca,'xticklabel',finalMatNames,'xticklabelRotation',45)
    %     ylabel('% of units')
    
    
%     CFMat = CFMat(logical(finalMat2(:,4)) & logical(finalMat2(:,6)) & logical(finalMat2(:,7)),:);
    
    % Plot each penetration
    % igure;
    % boxplot(CFMat(:,[1 4 6 7]))
    % ylabel('CF (Hz)');
    % set(gca,'xticklabel',finalMatNames([1 4 6 7]),'xticklabelRotation',45)
    % 
    % % Plot difference in prefered F0 in octave
    % idx = [4 6 7];
    % for i = 1:3,
    % %     d = log2(CFMat(:,1)./CFMat(:,idx(i)));
    %     figure;
    % %     hist(d,10)
    %     scatterQ(CFMat(:,1),CFMat(:,idx(i)),1);
    % end

    
end
cellProp = sum(finalMat2)./size(finalMat2,1)*100;
% Plot average over all implantations
figure;
bar(1:6,cellProp);
set(gca,'xticklabel',finalMatNames,'xticklabelRotation',45)
ylabel('% of units')
set(gca,'FontWeight','bold','YGrid','on')

% Corrected matrix,
% Remove cases where 1 1 0 0 1
finalMat3 = finalMat2;
for i = 2:6,
    idx = ~logical(finalMat3(:,i));
    finalMat3(idx,i+1:end) = 0;
end
cellProp = sum(finalMat3)./size(finalMat3,1)*100;
figure;
hold on;
bar(1,cellProp(1),'facecolor','k');
for i = 2:6,
bar(i,cellProp(i),'facecolor',col(i-1,:));
end
set(gca,'xtick',1:6,'xticklabel',finalMatNames,'xticklabelRotation',45)
ylabel('% of pitch units')
set(gca,'FontWeight','bold','YGrid','on')


% Plot difference in prefered F0 between all Harm and jitter 0
allCFMat2 = allCFMat(logical(finalMat3(:,2)),:);
figure;
scatterQ(allCFMat2(:,1),allCFMat2(:,2),1);
set(gca,'FontSize',20,'FontWeight','bold')
xlabel(['Pref. F0 - ' finalMatNames{1}]);
ylabel(['Pref. F0 - ' finalMatNames{2}]);
classical_stat_paired(allCFMat2(:,1),allCFMat2(:,2),1);


% Distribution of BF difference in octave between all harm and 0 jitter

d = abs(log2(allCFMat2(:,1)./allCFMat2(:,2)));
figure
boxplot(d)
set(gca,'FontSize',20,'FontWeight','bold','ygrid','on')
ylabel('\Delta F0 (octave)')
set(gca,'xtick','')
%
% for i = 1:6,
%     fprintf('%s: %2.2f%% Selective out of %d units\n',fList{i},allMat{i}.nU/allMat{i}.nAllU*100,allMat{i}.nAllU);
% end

%% Coding Strategy - Pooled all units
% Run Unit Selection first
%only figure 2 works but think that's how it is meant to be

stimNames1 = {'low' 'F0MaskLow'};
stimNames2 = {'allHarm' 'high' 'F0MaskHigh'};

a = 0.05;

finalMat = {};
CFMat = {};
finalMatNames = {'All Harm.' 'Resolved' 'Resolved Mskd' 'Unresolved' 'Unresolved Mskd' 'Both'};

h=figure;
col = lines;

for imp = 1:length(allDat),
    Dat = allDat{imp};
    uIdx = find(Dat.idxResp);
    CFValues = zeros(length(uIdx),2);
    nU = length(uIdx);
    nAllU = length(Dat.idxResp);
    
    finalMat = zeros(nU,5);
    CFMat = zeros(nU,5);
    
    for i = 1:nU,
        
        for j = 1:length(stimNames1),
            [meanMatLow(j,:), errMatLow(j,:), FListLow, pvalLow(j,:), CFlow(j)] = calc_TC(Dat.Data(uIdx(i)),gStim,gF0,stimNames1{j});
        end
        
        for j = 1:length(stimNames2),
            [meanMat(j,:), errMat(j,:), FList, pval(j,:), CF(j)] = calc_TC(Dat.Data(uIdx(i)),gStim,gF0,stimNames2{j});
        end
       
        finalMat(i,:) = [pval(1) pvalLow(1) pvalLow(2) pval(2) pval(3)] < a ;
        CFMat(i,:) = [CF(1) CFlow(1) CFlow(2) CF(2) CF(3)];
        
%         %&& logical(finalMat(i,4)) && logical(finalMat(i,5))
%         if logical(finalMat(i,2)) && ~logical(finalMat(i,3)) ,
%         figure(h); % Browsing for examples
%         clf
%         hold on;
%         for kk = 1,
% %             errorbar(FList,meanMat(kk,:),errMat(kk,:));
% % error_area( x,y,e,curve_col,area_col,linewidth,scale,X_ticks_label, Y_ticks_label)
%             error_area(FList,meanMat(kk,:),errMat(kk,:),col(kk,:),col(kk,:),1,'xlog',[200 400 800 1600 3200 6400]);
%         end
%         for kk = 1:2
%             error_area(FListLow,meanMatLow(kk,:),errMatLow(kk,:),col(kk+3,:),col(kk+3,:),1,'xlog',[200 400 800 1600 3200 6400]);
%         end
%         hold off;
% %         legend(finalMatNames(2:6));
% %         set(gca,'XScale','log')
% %         xlim([200 5000])
% %         set(gca,'XTick',[200 400 800 1600 3200])
%         xlabel('F0 (Hz)')
%         ylabel('Firing rate (spike/sec)')
%         set(gca,'FontSize',20,'FontWeight','bold')
%         pause();
%         end
%         
    end
    allMat{imp}.finalMat = finalMat;
    allMat{imp}.CFMat = CFMat;
    allMat{imp}.nU = nU;
    allMat{imp}.nAllU = nAllU;
end

allCFMat = [];
cellProp = [];
finalMat2 = [];
for imp = 1:length(allDat),
    
    idx = logical(allMat{imp}.finalMat(:,1));
%     idx = true(size(allMat{imp}.finalMat(:,1),1),1);
    finalMat2 = [finalMat2; allMat{imp}.finalMat(idx,:)];

    CFMat = allMat{imp}.CFMat(idx,:);
    allCFMat = [allCFMat; CFMat];
    
%     cellProp(imp,:) = sum(finalMat2)./size(finalMat2,1)*100
    
    % Plot each penetration
    %     figure;
    %     bar(sum(finalMat2)./size(finalMat2,1)*100)
    %     set(gca,'xticklabel',finalMatNames,'xticklabelRotation',45)
    %     ylabel('% of units')
    
    
%     CFMat = CFMat(logical(finalMat2(:,4)) & logical(finalMat2(:,6)) & logical(finalMat2(:,7)),:);
    
    % Plot each penetration
    % igure;
    % boxplot(CFMat(:,[1 4 6 7]))
    % ylabel('CF (Hz)');
    % set(gca,'xticklabel',finalMatNames([1 4 6 7]),'xticklabelRotation',45)
    % 
    % % Plot difference in prefered F0 in octave
    % idx = [4 6 7];
    % for i = 1:3,
    % %     d = log2(CFMat(:,1)./CFMat(:,idx(i)));
    %     figure;
    % %     hist(d,10)
    %     scatterQ(CFMat(:,1),CFMat(:,idx(i)),1);
    % end

    
end
finalMat2 = [finalMat2 finalMat2(:,3) & finalMat2(:,5)];
cellProp = sum(finalMat2)./size(finalMat2,1)*100;
% Plot average over all implantations
idxCol = [4 5 2 3];
figure;
hold on;
for i = 1:4,
    bar(i,cellProp(i+1),'facecolor',col(idxCol(i),:));
end
bar(5,cellProp(6),'facecolor','k');
hold off;
set(gca,'XTick',1:5,'xticklabel',finalMatNames(2:end),'xticklabelRotation',45)
ylabel('% of pitch units')
set(gca,'FontWeight','bold','YGrid','on')



% Get legend
figure;
plot(1,1,'color',col(1,:))
hold on
plot(1,1,'color',col(4,:))
plot(1,1,'color',col(5,:))

legend(finalMatNames([1 2 3]))
set(gca,'FontWeight' ,'bold')

%% Unit localisations
%works
finalDat = [];
for imp = 1:length(allDat),
    Dat = allDat{imp};
    uIdx = find(Dat.idxResp);
    CFValues = zeros(length(uIdx),2);
    nU = length(uIdx);
    nAllU = length(Dat.idxResp);
    nUTones = sum(arrayfun(@(x)(x.p(13)<0.05),Dat.DataAnova1));
    nUResp = sum(Dat.idxRespTtest);
    nUBoth = sum(arrayfun(@(x)(x.p(13)<0.05),Dat.DataAnova1) & Dat.idxResp);
    
    finalDat(imp,:) = [nAllU nUResp nU nUTones nUBoth];
    
end

PList = unique(LocPen);
finalDat2 = [];
finalDat2prc = [];
for i = 1:length(PList)
    idx = strcmp(LocPen,PList(i));
    finalDat2(i,:) = sum(finalDat(idx,:),1);
    finalDat2prc(i,:) = round(finalDat2(i,:) ./ finalDat2(i,1) *100,2);
end
disp(finalDat2)
disp(finalDat2prc)
disp(PList)

%% Find individual examples F0 tuning label('Pref. F0 (Hz)')
%works

stimNames1 = {'low', 'F0MaskLow'};
stimNames2 = {'high', 'F0MaskHigh', 'alt', 'rand', 'allHarm'};

col = lines;

a = 0.05;

finalMat = {};
CFMat = {};
finalMatNames = {'All Harm.' 'Resolved' 'Resolved Mskd' 'Unresolved' 'Unresolved Mskd' 'Unresolved alt' 'Unresolved rand'};

pairList = [8 13; 8 1; 8 11; 8 7; 8 10 ; 8 6;8 9; 8 12];
h = figure('position',[680   535   764   443]);

for imp = 1:length(allDat),
    Dat = allDat{imp};
    uIdx = find(Dat.idxResp);
    CFValues = zeros(length(uIdx),2);
    nU = length(uIdx);
    nAllU = length(Dat.idxResp);
    
    finalMat = zeros(nU,7);
    finalMatAH = zeros(nU,7);
    CFMat = zeros(nU,7);
    
    for i = 1:nU,
        
        for j = 1:length(stimNames1),
            [meanMatLow(j,:), errMatLow(j,:), FListLow, pvalLow(j,:), CFlow(j)] = calc_TC(Dat.Data(uIdx(i)),gStim,gF0,stimNames1{j});
        end
        
        for j = 1:length(stimNames2),
            [meanMat(j,:), errMat(j,:), FList, pval(j,:), CF(j)] = calc_TC(Dat.Data(uIdx(i)),gStim,gF0,stimNames2{j});
        end
        
        if sum(pval(5,:)<a)>0 %has tuning to all harm
            
            figure(h);
            error_area(FList,meanMat(5,:),errMat(5,:),col(1,:),col(1,:),1,'xlog',[200 400 800 1600 3200 6400]);
            hold on;
%             errorbar(FList,meanMat(1,:),errMat(1,:))
error_area(FList,meanMat(3,:),errMat(3,:),col(2,:),col(2,:),1,'xlog',[200 400 800 1600 3200 6400]);
error_area(FList,meanMat(4,:),errMat(4,:),col(3,:),col(3,:),1,'xlog',[200 400 800 1600 3200 6400]);
%             legend({'All Harm' 'Unresolved Harm' 'Alt Phase' 'Rand Phase'})
%             legend({'All Harm' 'Alt Phase' 'Rand Phase'})
%             set(gca,'xscale','log')
%             xlim([200 5000])
%             set(gca,'XTick',[200 400 800 1600 3200])
            xlabel('F0 (Hz)')
            ylabel('Firing rate (spike/sec)')
            set(gca,'FontSize',20,'FontWeight','bold')
            
            hold off;
            grid('on')
            fprintf('Impl : %d - Unit : %d\n',imp,i)
            pause();
            
            
            
            
            
            
            
            %         %&& logical(finalMat(i,4)) && logical(finalMat(i,5))
%         if logical(finalMat(i,2)) && ~logical(finalMat(i,3)) ,
%         figure(h); % Browsing for examples
%         clf
%         hold on;
%         for kk = 1,
% %             errorbar(FList,meanMat(kk,:),errMat(kk,:));
% % error_area( x,y,e,curve_col,area_col,linewidth,scale,X_ticks_label, Y_ticks_label)
%             error_area(FList,meanMat(kk,:),errMat(kk,:),col(kk,:),col(kk,:),1,'xlog',[200 400 800 1600 3200 6400]);
%         end
%         for kk = 1:2
%             error_area(FListLow,meanMatLow(kk,:),errMatLow(kk,:),col(kk+3,:),col(kk+3,:),1,'xlog',[200 400 800 1600 3200 6400]);
%         end
%         hold off;
% %         legend(finalMatNames(2:6));
% %         set(gca,'XScale','log')
% %         xlim([200 5000])
% %         set(gca,'XTick',[200 400 800 1600 3200])
%         xlabel('F0 (Hz)')
%         ylabel('Firing rate (spike/sec)')
%         set(gca,'FontSize',20,'FontWeight','bold')
%         pause();
%         end
            
        end
%                 
%         finalMat(i,:) = [(sum(pval(5,:)<a)>0) (sum(pvalLow(1,:)<a)>0) (sum(pvalLow(2,:)<a)>0) ...
%             (sum(pval(1,:)<a)>0) (sum(pval(2,:)<a)>0) (sum(pval(3,:)<a)>0) (sum(pval(4,:)<a)>0)];
%         CFMat(i,:) = [CF(5) CFlow(1:2) CF(1:4)];
    end
%     allMat{imp}.finalMat = finalMat;
%     allMat{imp}.finalMatAH = finalMatAH;
%     allMat{imp}.CFMat = CFMat;
%     allMat{imp}.nU = nU;
%     allMat{imp}.nAllU = nAllU;
end

%% Individual example Raster + F0 tuning
%works

stimNames1 = {'low', 'F0MaskLow'};
stimNames2 = {'high', 'F0MaskHigh', 'alt', 'rand', 'allHarm'};

stim = 'allHarm';
window = [0 0.5]; %seconds

a = 0.05;

finalMat = {};
CFMat = {};
finalMatNames = {'All Harm.' 'Resolved' 'Resolved Mskd' 'Unresolved' 'Unresolved Mskd' 'Unresolved alt' 'Unresolved rand'};

pairList = [8 13; 8 1; 8 11; 8 7; 8 10 ; 8 6;8 9; 8 12];
h1 = figure;
h2 = figure;

for imp = 1:length(allDat),
    changedstr = allDat{imp}.fName;
    changedstr = insertBefore(changedstr, 'Pitch', 'Good_'); %Dan added this to allow for new naming system
    load(['Spikes_' changedstr '.mat']); %if naming system changes back, change to allDat{imp}.fName
    
repList = unique(Y(:,5));

    Dat = allDat{imp};
    uIdx = find(Dat.idxResp);
    CFValues = zeros(length(uIdx),2);
    nU = length(uIdx);
    nAllU = length(Dat.idxResp);
    
    finalMat = zeros(nU,7);
    finalMatAH = zeros(nU,7);
    CFMat = zeros(nU,7);
    
    for i = 1:nU,
        
        for j = 1:length(stimNames1),
            [meanMatLow(j,:), errMatLow(j,:), FListLow, pvalLow(j,:), CFlow(j)] = calc_TC(Dat.Data(uIdx(i)),gStim,gF0,stimNames1{j});
        end
        
        for j = 1:length(stimNames2),
            [meanMat(j,:), errMat(j,:), FList, pval(j,:), CF(j)] = calc_TC(Dat.Data(uIdx(i)),gStim,gF0,stimNames2{j});
        end
        
        if sum(pval(5,:)<a)>0 %has tuning to all harm
            
            figure(h1);
            errorbar(FList,meanMat(5,:),errMat(5,:))
%             hold on;
%             errorbar(FList,meanMat(1,:),errMat(1,:))
%             errorbar(FList,meanMat(3,:),errMat(3,:))
%             errorbar(FList,meanMat(4,:),errMat(4,:))
%             legend({'All Harm' 'Unresolved Harm' 'Alt Phase' 'Rand Phase'})
%             legend({'All Harm' 'Alt Phase' 'Rand Phase'})
            set(gca,'xscale','log')
            xlim([200 5000])
            set(gca,'XTick',[200 400 800 1600 3200])
            xlabel('F0 (Hz)')
            ylabel('Firing rate (spike/sec)')
            set(gca,'FontSize',20,'FontWeight','bold')
            
%             hold off;


    Unit = allDat{imp}.Units(i);
    Nspikes = zeros(length(repList),length(Flist));
    SpikeTimes = {};
    c = 1;
    for f = 1:length(Flist),
        for r = 1:length(repList),
            stimNum = find(strcmp(type,stim) & (F0==Flist(f)));
            if isempty(stimNum),
                fprintf('Stim : %s , F0: %d',stim,Flist(f));
                continue;
            end
            idx = (Y(:,4) == stimNum) & (Y(:,3) == Unit) & (Y(:,5)==repList(r)) & (Y(:,2) > window(1)) & (Y(:,2) < window(2));
            Nspikes(r,f) = sum(idx);
            SpikeTimes{r,f} = Y(idx,2);
            c = c+1;
        end
    end
    Nspikes = Nspikes ./ diff(window);
    
    figure(h2);
%     clf
    rasterQ(SpikeTimes);
%     fprintf('Unit: %d - %s',Unit,stim)
    xlabel('Time (second)')
    ylabel('F0 (Hz)')
    set(gca,'YTickLabel',Flist)
    set(gca,'FontSize',16,'FontWeight','bold')

%     saveas(h2,sprintf('U%d_%s_Raster



            
            fprintf('Impl : %d - Unit : %d\n',imp,i)

            pause();
            
            
        end
        
        
    end
    
end



% 
% Unit = 330;
% stimList = unique(type);
% Flist = unique(F0);
% UnitList = unique(Y(:,3));
% repList = unique(Y(:,5));
% window = [0 0.5]; %seconds
% 
% h1 = figure;
% h2 = figure;
% 
% for i = 1:length(stimList),
%     stim = stimList{i};
%     Nspikes = zeros(length(repList),length(Flist));
%     SpikeTimes = {};
%     c = 1;
%     for f = 1:length(Flist),
%         for r = 1:length(repList),
%             stimNum = find(strcmp(type,stim) & (F0==Flist(f)));
%             if isempty(stimNum),
%                 fprintf('Stim : %s , F0: %d',stim,Flist(f));
%                 continue;
%             end
%             idx = (Y(:,4) == stimNum) & (Y(:,3) == Unit) & (Y(:,5)==repList(r)) & (Y(:,2) > window(1)) & (Y(:,2) < window(2));
%             Nspikes(r,f) = sum(idx);
%             SpikeTimes{c} = Y(idx,2);
%             c = c+1;
%         end
%     end
%     Nspikes = Nspikes ./ diff(window);
% %     figure(h1);
% %     clf
% %     errorbar(Flist,mean(Nspikes),ste(Nspikes))
% %     title(sprintf('Unit: %d - %s',Unit,stim))
% %     set(gca,'XScale','log')
% %     xlabel('F0 (Hz)')
% %     ylabel('Evoked Firing Rate (spike/sec)');
% %     saveas(h1,sprintf('U%d_%s_Selectivity.png',Unit,stim));
%     
%     figure(h2);
%     clf
%     hold on
%     for k = 1:length(SpikeTimes)
%         scatter(SpikeTimes{k},ones(length(SpikeTimes{k}),1) * k,'marker','x','markerfacecolor','k','markeredgecolor','k')
%     end
%     y = length(repList):length(repList):length(repList)*(length(Flist));
%     line(repmat(window,length(y),1)',[y' y']','color','r');
%     hold off
%     title(sprintf('Unit: %d - %s',Unit,stim))
%     xlabel('Time relative to stimulus onset (second)')
%     ylabel('F0 (Hz)')
%     set(gca,'YTick',y)
%     set(gca,'YTickLabel',Flist)
% %     saveas(h2,sprintf('U%d_%s_Raster.png',Unit,stim));
%     
% %     pause();
% end

%% Coding Strategy - Teaser version
% Run Unit Selection first

stimNames1 = {'low' 'F0MaskLow'};
stimNames2 = {'allHarm' 'high' 'F0MaskHigh'};

a = 0.05;

finalMat = {};
CFMat = {};
finalMatNames = {'All Harm.' 'Resolved' 'Resolved Mskd' 'Unresolved' 'Unresolved Mskd' 'Both'};

h=figure;
col = lines;

for imp = 1:length(allDat),
    Dat = allDat{imp};
    uIdx = find(Dat.idxResp);
    CFValues = zeros(length(uIdx),2);
    nU = length(uIdx);
    nAllU = length(Dat.idxResp);
    
    finalMat = zeros(nU,5);
    CFMat = zeros(nU,5);
    
    for i = 1:nU,
        
        for j = 1:length(stimNames1),
            [meanMatLow(j,:), errMatLow(j,:), FListLow, pvalLow(j,:), CFlow(j)] = calc_TC(Dat.Data(uIdx(i)),gStim,gF0,stimNames1{j});
        end
        
        for j = 1:length(stimNames2),
            [meanMat(j,:), errMat(j,:), FList, pval(j,:), CF(j)] = calc_TC(Dat.Data(uIdx(i)),gStim,gF0,stimNames2{j});
        end
       
        finalMat(i,:) = [pval(1) pvalLow(1) pvalLow(2) pval(2) pval(3)] < a ;
        CFMat(i,:) = [CF(1) CFlow(1) CFlow(2) CF(2) CF(3)];
        
        %&& logical(finalMat(i,4)) && logical(finalMat(i,5))
        % logical(finalMat(i,2)) && ~logical(finalMat(i,3))
        if true,
        figure(h); % Browsing for examples
        clf
        hold on;
        for kk = [1 3],
%             errorbar(FList,meanMat(kk,:),errMat(kk,:));
% error_area( x,y,e,curve_col,area_col,linewidth,scale,X_ticks_label, Y_ticks_label)
            error_area(FList,meanMat(kk,:),errMat(kk,:),col(kk,:),col(kk,:),1,'xlog',[200 400 800 1600 3200 6400]);
        end
        for kk = 2
            error_area(FListLow,meanMatLow(kk,:),errMatLow(kk,:),col(kk+3,:),col(kk+3,:),1,'xlog',[200 400 800 1600 3200 6400]);
        end
        hold off;
%         legend(finalMatNames(2:6));
%         set(gca,'XScale','log')
%         xlim([200 5000])
%         set(gca,'XTick',[200 400 800 1600 3200])
        xlabel('F0 (Hz)')
        ylabel('Firing rate (spike/sec)')
        set(gca,'FontSize',20,'FontWeight','bold')
        pause();
        end
%         
    end
    allMat{imp}.finalMat = finalMat;
    allMat{imp}.CFMat = CFMat;
    allMat{imp}.nU = nU;
    allMat{imp}.nAllU = nAllU;
end

figure;
plot(1,1,'color',col(1,:))
hold on;
plot(1,1,'color',col(3,:))
plot(1,1,'color',col(5,:));
legend(finalMatNames([1 4 2]))
set(gca,'FontWeight','bold')
