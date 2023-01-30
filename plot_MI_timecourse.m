
%% Plot the MI timecourse for each stimulus type both for each neuron and averaged for each penetration
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, January 2023

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

load('unit_MIs_30ms.mat')

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