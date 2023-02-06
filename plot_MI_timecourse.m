
%% Plot the MI timecourse for each stimulus type both for each neuron and averaged for each penetration
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, January 2023

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

% load('unit_MIs_30ms.mat')

% for each penetration
for pen = 1:length(unit_MIs)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit_MIs{pen,1} '/tmp/Spikes_' unit_MIs{pen,1} '_' unit_MIs{pen,2} '_Good_Pitch.mat']);

    allUnit_MI_tcs = unit_MIs{pen,3};

    units = unique(Y(:,3));

    all_lows = [];
    all_highs = [];

    for uu = 1:length(units)

        if ~isempty(find(allUnit_MI_tcs(uu,11,:), 1))
            all_lows = [all_lows; allUnit_MI_tcs(uu,11,:)];
        end

        if ~isempty(find(allUnit_MI_tcs(uu,10,:),1))
            all_highs = [all_highs; allUnit_MI_tcs(uu,10,:)];
        end


%         if ~isempty(find(allUnit_MI_tcs(uu,[1 11 10 13],:),1))
% 
%         
%             % plot CT0, low, high, tone, allHarm
%             %       1    11   10    13      8
%             figure('Position',[1900 500 1800 1200])
%     
%             %%%%%% plot CT0 %%%%%%
%             subplot(5,1,1)
%             plot(squeeze(allUnit_MI_tcs(uu,1,:)),'k')
%     
%             subplot(5,1,2)
%             plot(squeeze(allUnit_MI_tcs(uu,11,:)),'r')
%     
%             subplot(5,1,3)
%             plot(squeeze(allUnit_MI_tcs(uu,10,:)),'b')
%     
%             subplot(5,1,4)
%             plot(squeeze(allUnit_MI_tcs(uu,13,:)),'g')
%     
%             subplot(5,1,5)
%             plot(squeeze(allUnit_MI_tcs(uu,8,:)),'c')
%     
%             pause
%         end

    end % ends unit loop

    figure('Position',[1900 500 1800 1200])
    avg_lows = squeeze(mean(all_lows,'omitnan'));
    avg_highs = squeeze(mean(all_highs,'omitnan'));
    plot(avg_lows,'linewidth',3); hold on
    plot(avg_highs,'linewidth',3)

    axis tight;
    xticks(0:5:30)
    xticklabels(15:15*5:500)

    legend('Low','High')

    a=1;

end % ends penetration loop



%% Compare the different time bins

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13


load('unit_MIs_50ms.mat')
uMI50ms = unit_MIs;

load('unit_MIs_30ms.mat')
uMI30ms = unit_MIs;

load('unit_MIs_20ms.mat')
uMI20ms = unit_MIs;

load('unit_MIs_50ms_fc.mat')
uMI50msfc = unit_MIs;

load('unit_MIs_30ms_fc.mat')
uMI30msfc = unit_MIs;

load('unit_MIs_20ms_fc.mat')
uMI20msfc = unit_MIs;

clear unit_MIs;

for pen = 1:length(uMI50ms)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' uMI50ms{pen,1} '/tmp/Spikes_' uMI50ms{pen,1} '_' uMI50ms{pen,2} '_Good_Pitch.mat']);
    
    MI50 = uMI50ms{pen,3};
    MI30 = uMI30ms{pen,3};
    MI20 = uMI20ms{pen,3};
    MI50fc = uMI50msfc{pen,3};
    MI30fc = uMI30msfc{pen,3};
    MI20fc = uMI20msfc{pen,3};

    low_not_sensitive = find(sensitivity(:,11)==0);

    low50 = squeeze(MI50(:,11,:));
    low50(low_not_sensitive,:) = [];
    avg_low50 = mean(low50);
    
    low30 = squeeze(MI30(:,11,:));
    low30(low_not_sensitive,:) = [];
    avg_low30 = mean(low30);

    low20 = squeeze(MI20(:,11,:));
    low20(low_not_sensitive,:) = [];
    avg_low20 = mean(low20);

    low50fc = squeeze(MI50fc(:,11,:));
    low50fc(low_not_sensitive,:) = [];
    avg_low50fc = mean(low50fc);

    low30fc = squeeze(MI30fc(:,11,:));
    low30fc(low_not_sensitive,:) = [];
    avg_low30fc = mean(low30fc);

    low20fc = squeeze(MI20fc(:,11,:));
    low20fc(low_not_sensitive,:) = [];
    avg_low20fc = mean(low20fc);
    

    %%%%% PLOT %%%%%

    figure('Position',[1900 500 1800 1200])
    hold on
    plot(1:25:475,avg_low50,'k','linewidth',3)
    plot(1:15:475,avg_low30,'b','linewidth',3)
    plot(1:10:485,avg_low20,'r','linewidth',3)

    plot(1:25:475,avg_low50fc,'k--','linewidth',4)
    plot(1:15:475,avg_low30fc,'b--','linewidth',4)
    plot(1:10:485,avg_low20fc,'r--','linewidth',4)

    legend('50','30','20')

    xticks(0:50:500)
    axis tight

    set(gca,'fontsize',22)

    pause;
end % ends recording loop


%% Plot MI timecourses by region

lowA1 = [11 15 16 17];
lowA1_N = [1 2 4];

highA1 = [12 14];
highA1_N = 3;

lowAAF = [9 10 13 18 20];
lowAAF_N = 6;

highAAF = 19;
highAAF_N = 8;

PPF = [5 7];

load('unit_MIs_30ms.mat')

stims = {'CT0','CT40','allHarm','tone','low','high'};
stimIDX = [1 4 8 13 11 10];


%%%%%% low A1 plot %%%%%

figure('Position',[1900 500 1800 1200])
sgtitle('Low A1','fontsize',32)

for ss = 1:length(stims)

    all = [];
    for pen = lowA1

        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit_MIs{pen,1} '/tmp/Spikes_' unit_MIs{pen,1} '_' unit_MIs{pen,2} '_Good_Pitch.mat']);

        MIs = unit_MIs{pen,3};
    
        not_sensitive = find(sensitivity(:,stimIDX(ss))==0);
    
        stim_MIs = squeeze(MIs(:,stimIDX(ss),:));
        stim_MIs(not_sensitive,:) = [];
        stim_MIs(:,end-6:end) = [];
        all = [all; stim_MIs];
    end

    for pen = lowA1_N

        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit_MIs{pen,1} '/tmp/Spikes_' unit_MIs{pen,1} '_' unit_MIs{pen,2} '_Good_Pitch.mat']);

        MIs = unit_MIs{pen,3};
    
        not_sensitive = find(sensitivity(:,stimIDX(ss))==0);
    
        stim_MIs = squeeze(MIs(:,stimIDX(ss),:));
        stim_MIs(not_sensitive,:) = [];
        stim_MIs(:,14:20) = [];
        all = [all; stim_MIs];
    end

    avg_MI = mean(all);
    std_MI = std(all);

    subplot(3,2,ss)
    shadedErrorBar(0:15:370,avg_MI,std_MI,{'b','linewidth',3})
    hold on
    ymax = 0.18;
    p = patch([0 0 200 200],[0 ymax ymax 0],'k','FaceAlpha',0.1,'edgecolor','none');
    axis tight
    ylim([0 ymax])
    yticks(0:0.06:ymax)

    if ss == length(stims)-1
        xlabel('ms since stimulus onset')
        ylabel('MI')
    end

    title(stims(ss))

    set(gca,'fontsize',20)

end


%%%%%% low AAF plot %%%%%

figure('Position',[1900 500 1800 1200])
sgtitle('Low AAF','fontsize',32)

for ss = 1:length(stims)

    all = [];
    for pen = lowAAF

        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit_MIs{pen,1} '/tmp/Spikes_' unit_MIs{pen,1} '_' unit_MIs{pen,2} '_Good_Pitch.mat']);

        MIs = unit_MIs{pen,3};
    
        not_sensitive = find(sensitivity(:,stimIDX(ss))==0);
    
        stim_MIs = squeeze(MIs(:,stimIDX(ss),:));
        stim_MIs(not_sensitive,:) = [];
        stim_MIs(:,end-6:end) = [];
        all = [all; stim_MIs];
    end

    for pen = lowAAF_N

        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit_MIs{pen,1} '/tmp/Spikes_' unit_MIs{pen,1} '_' unit_MIs{pen,2} '_Good_Pitch.mat']);

        MIs = unit_MIs{pen,3};
    
        not_sensitive = find(sensitivity(:,stimIDX(ss))==0);
    
        stim_MIs = squeeze(MIs(:,stimIDX(ss),:));
        stim_MIs(not_sensitive,:) = [];
        stim_MIs(:,14:20) = [];
        all = [all; stim_MIs];
    end

    avg_MI = mean(all);
    std_MI = std(all);

    subplot(3,2,ss)
    shadedErrorBar(0:15:370,avg_MI,std_MI,{'b','linewidth',3})
    hold on
    ymax = 0.33;
    p = patch([0 0 200 200],[0 ymax ymax 0],'k','FaceAlpha',0.1,'edgecolor','none');
    axis tight
    ylim([0 ymax])
    yticks(0:.1:ymax)

    if ss == length(stims)-1
        xlabel('ms since stimulus onset')
        ylabel('MI')
    end
    
    title(stims(ss))

    set(gca,'fontsize',20)

end



%%%%%% high A1 plot %%%%%

figure('Position',[1900 500 1800 1200])
sgtitle('High A1','fontsize',32)

for ss = 1:length(stims)

    all = [];
    for pen = highA1

        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit_MIs{pen,1} '/tmp/Spikes_' unit_MIs{pen,1} '_' unit_MIs{pen,2} '_Good_Pitch.mat']);

        MIs = unit_MIs{pen,3};
    
        not_sensitive = find(sensitivity(:,stimIDX(ss))==0);
    
        stim_MIs = squeeze(MIs(:,stimIDX(ss),:));
        stim_MIs(not_sensitive,:) = [];
        stim_MIs(:,end-6:end) = [];
        all = [all; stim_MIs];
    end

    for pen = highA1_N

        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit_MIs{pen,1} '/tmp/Spikes_' unit_MIs{pen,1} '_' unit_MIs{pen,2} '_Good_Pitch.mat']);

        MIs = unit_MIs{pen,3};
    
        not_sensitive = find(sensitivity(:,stimIDX(ss))==0);
    
        stim_MIs = squeeze(MIs(:,stimIDX(ss),:));
        stim_MIs(not_sensitive,:) = [];
        stim_MIs(:,14:20) = [];
        all = [all; stim_MIs];
    end

    avg_MI = mean(all);
    std_MI = std(all);

    subplot(3,2,ss)
    shadedErrorBar(0:15:370,avg_MI,std_MI,{'b','linewidth',3})
    hold on
    ymax = 0.23;
    p = patch([0 0 200 200],[0 ymax ymax 0],'k','FaceAlpha',0.1,'edgecolor','none');
    axis tight
    ylim([0 ymax])
    yticks(0:0.08:ymax)

    if ss == length(stims)-1
        xlabel('ms since stimulus onset')
        ylabel('MI')
    end
    
    title(stims(ss))

    set(gca,'fontsize',20)

end


%%%%%% high AAF plot %%%%%

figure('Position',[1900 500 1800 1200])
sgtitle('High AAF','fontsize',32)

for ss = 1:length(stims)

    all = [];
    for pen = highAAF

        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit_MIs{pen,1} '/tmp/Spikes_' unit_MIs{pen,1} '_' unit_MIs{pen,2} '_Good_Pitch.mat']);

        MIs = unit_MIs{pen,3};
    
        not_sensitive = find(sensitivity(:,stimIDX(ss))==0);
    
        stim_MIs = squeeze(MIs(:,stimIDX(ss),:));
        stim_MIs(not_sensitive,:) = [];
        stim_MIs(:,end-6:end) = [];
        all = [all; stim_MIs];
    end

    for pen = highAAF_N

        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit_MIs{pen,1} '/tmp/Spikes_' unit_MIs{pen,1} '_' unit_MIs{pen,2} '_Good_Pitch.mat']);

        MIs = unit_MIs{pen,3};
    
        not_sensitive = find(sensitivity(:,stimIDX(ss))==0);
    
        stim_MIs = squeeze(MIs(:,stimIDX(ss),:));
        stim_MIs(not_sensitive,:) = [];
        stim_MIs(:,14:20) = [];
        all = [all; stim_MIs];
    end

    avg_MI = mean(all);
    std_MI = std(all);

    subplot(3,2,ss)
    shadedErrorBar(0:15:370,avg_MI,std_MI,{'b','linewidth',3})
    hold on
    ymax = 0.32;
    p = patch([0 0 200 200],[0 ymax ymax 0],'k','FaceAlpha',0.1,'edgecolor','none');
    axis tight
    ylim([0 ymax])
    yticks(0:.1:ymax)

    if ss == length(stims)-1
        xlabel('ms since stimulus onset')
        ylabel('MI')
    end
    
    title(stims(ss))

    set(gca,'fontsize',20)

end


%%%%%% PPF plot %%%%%

figure('Position',[1900 500 1800 1200])
sgtitle('PPF','fontsize',32)

for ss = 1:length(stims)

    all = [];

    for pen = PPF

        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit_MIs{pen,1} '/tmp/Spikes_' unit_MIs{pen,1} '_' unit_MIs{pen,2} '_Good_Pitch.mat']);

        MIs = unit_MIs{pen,3};
    
        not_sensitive = find(sensitivity(:,stimIDX(ss))==0);
    
        stim_MIs = squeeze(MIs(:,stimIDX(ss),:));
        stim_MIs(not_sensitive,:) = [];
        stim_MIs(:,14:20) = [];
        all = [all; stim_MIs];
    end

    avg_MI = mean(all,'omitnan');
    std_MI = std(all,'omitnan');

    if length(std_MI)==1
        avg_MI = all;
        std_MI = zeros(length(avg_MI),1);
    end

    subplot(3,2,ss)
    shadedErrorBar(0:15:370,avg_MI,std_MI,{'b','linewidth',3})
    hold on
    ymax = 0.22;
    p = patch([0 0 200 200],[0 ymax ymax 0],'k','FaceAlpha',0.1,'edgecolor','none');
    axis tight
    ylim([0 ymax])
    yticks(0:0.08:ymax)

    if ss == length(stims)-1
        xlabel('ms since stimulus onset')
        ylabel('MI')
    end
    
    title(stims(ss))

    set(gca,'fontsize',20)

end



%% Plot the timings of first and second peaks


figure('Position',[1900 500 1800 1200])

lowA1 = [11 15 16 17];
lowA1_N = [1 2 4];

highA1 = [12 14];
highA1_N = 3;

lowAAF = [9 10 13 18 20];
lowAAF_N = 6;

highAAF = 19;
highAAF_N = 8;

PPF = [5 7];

load('unit_MIs_30ms.mat')

stims = {'CT0','CT40','allHarm','tone','low','high'};
stimIDX = [1 4 8 13 11 10];
colors = colormap(hsv(length(stims)));

timebins = 0:15:500;
peak_thresh = 0.8;

peak_times = zeros(length(stims),2);
err = zeros(length(stims),2);

% highAAF = 9:20;
% highAAF_N = 1:8;

all_peak_times = cell(length(stims),1);

%%%%%%%%%%%%%% high AAF %%%%%%%%%%%%%%%%%%%%%%
for ss = 1:length(stims)

    tp1 = [];
    tp2 = [];

    for pen = highAAF

        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit_MIs{pen,1} '/tmp/Spikes_' unit_MIs{pen,1} '_' unit_MIs{pen,2} '_Good_Pitch.mat']);

        MIs = unit_MIs{pen,3};
    
        not_sensitive = find(sensitivity(:,stimIDX(ss))==0);
    
        stim_MIs = squeeze(MIs(:,stimIDX(ss),:));
        stim_MIs(not_sensitive,:) = [];
        stim_MIs(:,end-6:end) = [];

        for uu = 1:size(stim_MIs,1)
            [nPeaks, peakIDX] = count_peaks(stim_MIs(uu,:),peak_thresh);
            if nPeaks>1
                tp1 = [tp1; timebins(peakIDX(1))];
                tp2 = [tp2; timebins(peakIDX(2))-200];
            else
                tp1 = [tp1; timebins(peakIDX(1))];
            end
        end
    end

    
    for pen = highAAF_N

        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit_MIs{pen,1} '/tmp/Spikes_' unit_MIs{pen,1} '_' unit_MIs{pen,2} '_Good_Pitch.mat']);

        MIs = unit_MIs{pen,3};
    
        not_sensitive = find(sensitivity(:,stimIDX(ss))==0);
    
        stim_MIs = squeeze(MIs(:,stimIDX(ss),:));
        stim_MIs(not_sensitive,:) = [];
        stim_MIs(:,14:20) = [];
%         stim_MIs(:,21:31) = [];

        for uu = 1:size(stim_MIs,1)
            [nPeaks, peakIDX] = count_peaks(stim_MIs(uu,:),peak_thresh);

            if nPeaks > 1
                tp1 = [tp1; timebins(peakIDX(1))];
                tp2 = [tp2; timebins(peakIDX(2))-200];
            else
                tp1 = [tp1; timebins(peakIDX(1))];
            end
        end
        
    end

    all_peak_times{ss} = tp1;

    peak_times(ss,1) = mean(tp1);
    peak_times(ss,2) = mean(tp2);

    err(ss,1) = std(tp1);
    err(ss,2) = std(tp2);

end

peak_times = peak_times';
err = err';

subplot(2,2,1)
b = bar(peak_times,'Facecolor','flat');
for k = 1:length(stims)
    b(k).CData = colors(k,:);
end
hold on;
ngroups = size(peak_times,1);
nbars = size(peak_times,2);
groupwidth = min(0.8, nbars/(nbars+1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth/(2*nbars);
    errorbar(x,peak_times(:,i),err(:,i),'k.');
end

title('High AAF')

% ylim([0 360])
yticks([0:100:360])
set(gca,'Fontsize',24)
xlabel('Peak')


%%%%%%%%%%%%%% high A1 %%%%%%%%%%%%%%%%%%%%%%
peak_times = zeros(length(stims),2);
err = zeros(length(stims),2);

for ss = 1:length(stims)

    tp1 = [];
    tp2 = [];

    for pen = highA1

        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit_MIs{pen,1} '/tmp/Spikes_' unit_MIs{pen,1} '_' unit_MIs{pen,2} '_Good_Pitch.mat']);

        MIs = unit_MIs{pen,3};
    
        not_sensitive = find(sensitivity(:,stimIDX(ss))==0);
    
        stim_MIs = squeeze(MIs(:,stimIDX(ss),:));
        stim_MIs(not_sensitive,:) = [];
        stim_MIs(:,end-6:end) = [];

        for uu = 1:size(stim_MIs,1)
            [nPeaks, peakIDX] = count_peaks(stim_MIs(uu,:),peak_thresh);
%             clf;
%             plot(stim_MIs(uu,:));
%             hold on
%             scatter(peakIDX,stim_MIs(uu,peakIDX),100)
%             peakIDX
%             a=1;

            if nPeaks>1
                tp1 = [tp1; timebins(peakIDX(1))];
                tp2 = [tp2; timebins(peakIDX(2))-200];
            else
                tp1 = [tp1; timebins(peakIDX(1))];
            end
        end

    end



    for pen = highA1_N

        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit_MIs{pen,1} '/tmp/Spikes_' unit_MIs{pen,1} '_' unit_MIs{pen,2} '_Good_Pitch.mat']);

        MIs = unit_MIs{pen,3};
    
        not_sensitive = find(sensitivity(:,stimIDX(ss))==0);
    
        stim_MIs = squeeze(MIs(:,stimIDX(ss),:));
        stim_MIs(not_sensitive,:) = [];
        stim_MIs(:,14:20) = [];

        for uu = 1:size(stim_MIs,1)
            [nPeaks, peakIDX] = count_peaks(stim_MIs(uu,:),peak_thresh);

            if nPeaks>1
                tp1 = [tp1; timebins(peakIDX(1))];
                tp2 = [tp2; timebins(peakIDX(2))-200];
            else
                tp1 = [tp1; timebins(peakIDX(1))];
            end
        end    
    end

    peak_times(ss,1) = mean(tp1);
    peak_times(ss,2) = mean(tp2);

    err(ss,1) = std(tp1);
    err(ss,2) = std(tp2);

end

peak_times = peak_times';
err = err';

subplot(2,2,2)
b = bar(peak_times,'Facecolor','flat');
for k = 1:length(stims)
    b(k).CData = colors(k,:);
end
hold on;
ngroups = size(peak_times,1);
nbars = size(peak_times,2);
groupwidth = min(0.8, nbars/(nbars+1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth/(2*nbars);
    errorbar(x,peak_times(:,i),err(:,i),'k.');
end

title('High A1')

% ylim([0 360])
yticks([0:100:360])
set(gca,'Fontsize',24)
xlabel('Peak')



%%%%%%%%%%%%%% low AAF %%%%%%%%%%%%%%%%%%%%%%

peak_times = zeros(length(stims),2);
err = zeros(length(stims),2);

for ss = 1:length(stims)

    tp1 = [];
    tp2 = [];

    for pen = lowAAF

        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit_MIs{pen,1} '/tmp/Spikes_' unit_MIs{pen,1} '_' unit_MIs{pen,2} '_Good_Pitch.mat']);

        MIs = unit_MIs{pen,3};
    
        not_sensitive = find(sensitivity(:,stimIDX(ss))==0);
    
        stim_MIs = squeeze(MIs(:,stimIDX(ss),:));
        stim_MIs(not_sensitive,:) = [];
        stim_MIs(:,end-6:end) = [];

        for uu = 1:size(stim_MIs,1)
            [nPeaks, peakIDX] = count_peaks(stim_MIs(uu,:),peak_thresh);
%             clf;
%             plot(stim_MIs(uu,:));
%             hold on
%             scatter(peakIDX,stim_MIs(uu,peakIDX),100)
%             peakIDX
%             a=1;

            if nPeaks>1
                tp1 = [tp1; timebins(peakIDX(1))];
                tp2 = [tp2; timebins(peakIDX(2))-200];
            else
                tp1 = [tp1; timebins(peakIDX(1))];
            end
        end
    end



    for pen = lowAAF_N

        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit_MIs{pen,1} '/tmp/Spikes_' unit_MIs{pen,1} '_' unit_MIs{pen,2} '_Good_Pitch.mat']);

        MIs = unit_MIs{pen,3};
    
        not_sensitive = find(sensitivity(:,stimIDX(ss))==0);
    
        stim_MIs = squeeze(MIs(:,stimIDX(ss),:));
        stim_MIs(not_sensitive,:) = [];
        stim_MIs(:,14:20) = [];

        for uu = 1:size(stim_MIs,1)
            [nPeaks, peakIDX] = count_peaks(stim_MIs(uu,:),peak_thresh);

            if nPeaks>1
                tp1 = [tp1; timebins(peakIDX(1))];
                tp2 = [tp2; timebins(peakIDX(2))-200];
            else
                tp1 = [tp1; timebins(peakIDX(1))];
            end
        end


        
    end

    peak_times(ss,1) = mean(tp1);
    peak_times(ss,2) = mean(tp2);

    err(ss,1) = std(tp1);
    err(ss,2) = std(tp2);

end

peak_times = peak_times';
err = err';

subplot(2,2,3)
b = bar(peak_times,'Facecolor','flat');
for k = 1:length(stims)
    b(k).CData = colors(k,:);
end
hold on;
ngroups = size(peak_times,1);
nbars = size(peak_times,2);
groupwidth = min(0.8, nbars/(nbars+1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth/(2*nbars);
    errorbar(x,peak_times(:,i),err(:,i),'k.');
end

title('Low AAF')

% ylim([0 360])
yticks([0:100:360])
set(gca,'Fontsize',24)
xlabel('Peak')
ylabel('Timing relative to stim onset')



%%%%%%%%%%%%%% low A1 %%%%%%%%%%%%%%%%%%%%%%

peak_times = zeros(length(stims),2);
err = zeros(length(stims),2);

for ss = 1:length(stims)

    tp1 = [];
    tp2 = [];

    for pen = lowA1

        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit_MIs{pen,1} '/tmp/Spikes_' unit_MIs{pen,1} '_' unit_MIs{pen,2} '_Good_Pitch.mat']);

        MIs = unit_MIs{pen,3};
    
        not_sensitive = find(sensitivity(:,stimIDX(ss))==0);
    
        stim_MIs = squeeze(MIs(:,stimIDX(ss),:));
        stim_MIs(not_sensitive,:) = [];
        stim_MIs(:,end-6:end) = [];

        for uu = 1:size(stim_MIs,1)
            [nPeaks, peakIDX] = count_peaks(stim_MIs(uu,:),peak_thresh);
%             clf;
%             plot(stim_MIs(uu,:));
%             hold on
%             scatter(peakIDX,stim_MIs(uu,peakIDX),100)
%             peakIDX
%             a=1;

            if nPeaks>1
                tp1 = [tp1; timebins(peakIDX(1))];
                tp2 = [tp2; timebins(peakIDX(2))-200];
            else
                tp1 = [tp1; timebins(peakIDX(1))];
            end

        end

    end


    for pen = lowA1_N

        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit_MIs{pen,1} '/tmp/Spikes_' unit_MIs{pen,1} '_' unit_MIs{pen,2} '_Good_Pitch.mat']);

        MIs = unit_MIs{pen,3};
    
        not_sensitive = find(sensitivity(:,stimIDX(ss))==0);
    
        stim_MIs = squeeze(MIs(:,stimIDX(ss),:));
        stim_MIs(not_sensitive,:) = [];
        stim_MIs(:,14:20) = [];

        for uu = 1:size(stim_MIs,1)
            [nPeaks, peakIDX] = count_peaks(stim_MIs(uu,:),peak_thresh);

            if nPeaks>1
                tp1 = [tp1; timebins(peakIDX(1))];
                tp2 = [tp2; timebins(peakIDX(2))-200];
            else
                tp1 = [tp1; timebins(peakIDX(1))];
            end
        
        end

    end

    peak_times(ss,1) = mean(tp1);
    peak_times(ss,2) = mean(tp2);

    err(ss,1) = std(tp1);
    err(ss,2) = std(tp2);

end

peak_times = peak_times';
err = err';

subplot(2,2,4)
b = bar(peak_times,'Facecolor','flat');
for k = 1:length(stims)
    b(k).CData = colors(k,:);
end
hold on;
ngroups = size(peak_times,1);
nbars = size(peak_times,2);
groupwidth = min(0.8, nbars/(nbars+1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth/(2*nbars);
    errorbar(x,peak_times(:,i),err(:,i),'k.');
end

% ylim([0 360])
yticks([0:100:360])
set(gca,'Fontsize',24)
xlabel('Peak')


title('Low A1')
legend(stims,'location','northwest')


%% Overlay all MI timecourses

load('unit_MIs_20ms.mat')

figure('Position',[1900 500 1800 1200])

% stims = [1 10 11 13];
% stim_names = {'CT0','high','low','tone'};

stim_names = {'CT0','CT40','allHarm','tone','low','high'};
stims = [1 4 8 13 11 10];


for pen = 1:length(unit_MIs)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit_MIs{pen,1} '/tmp/Spikes_' unit_MIs{pen,1} '_' unit_MIs{pen,2} '_Good_Pitch.mat']);

    allUnit_MI_tcs = unit_MIs{pen,3};

    units = unique(Y(:,3));

    for uu = 1:length(units)
        unit = units(uu);

        for ss = 1:length(stims)
            stim = stims(ss);

            if sensitivity(uu,stim)~=0
                subplot(3,2,ss)
                hold on
                if strcmp(unit_MIs{pen,1},'Noah')
                    MI_tcs = squeeze(allUnit_MI_tcs(uu,stim,:));
                    MI_tcs(21:31) = [];
                    plot(MI_tcs,'k')
                else
                    plot(squeeze(allUnit_MI_tcs(uu,stim,:)),'k')
                end
                axis tight
                
            end
        end

    end
end

for ss = 1:length(stims)
    subplot(3,2,ss)
    title(stim_names{ss})
    ylim([0 0.7])
    set(gca,'fontsize',22)
%     xticks(0:5:30)
%     xticklabels(0:75:500)
end