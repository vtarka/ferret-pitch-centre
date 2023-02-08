
%% Plots xxxxx
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, February 2023

load('all_unit_MIs_30ms_2.mat')

all_peaks = zeros(1311,13,2);
unit = 1;

for pen = 1:length(unit_MIs)

    MIs = unit_MIs{pen,3};
    active = unit_MIs{pen,4};

    [r,c] = find(~active);

    MIs(r,c,:) = -1;

    for uu = 1:size(MIs,1)
        
        for ss = 1:13
            
            if MIs(uu,ss,1) ~= -1
                [nPeaks, peakIDX] = count_peaks(MIs(uu,ss,:),0.8);

                if peakIDX(1) < 9
                    all_peaks(unit,ss,1) = 9 - peakIDX(1);
                end
                
                if strcmp(unit_MIs{pen,1},'Noah')
                    if ~isempty(find(peakIDX > 21, 1))
                        temp = peakIDX(peakIDX>21);
                        if temp(1) < 28
                            all_peaks(unit,ss,2) = 28 - temp(1);
                        end
                    end
                else
                    if ~isempty(find(peakIDX > 13, 1))
                        temp = peakIDX(peakIDX>13);
                        if temp(1) < 21
                            all_peaks(unit,ss,2) = 21 - temp(1);
                        end
                    end
                end
            end

        end % ends stim loop

        unit = unit + 1;

    end % ends unit loop
end % ends recording loop


%% plot

CT0_peaks = squeeze(all_peaks(:,4,:));
sortedCT0peaks = sortrows(CT0_peaks,[-1 -2]);

figure; imagesc(sortedCT0peaks)

xticks(1:2)
set(gca,'fontsize',24)
c=colorbar;
c.Ticks = linspace(0,8,8);
c.TickLabels = {'No Peak','105','90','75','60','45','30','15'};

xlabel('Peak after onset (1) or offset (2)')
ylabel('Neuron')
title('CT40')
ylabel(c,'ms after stim onset (1) or offset (2)','fontsize',24,'rotation',270)



%% plot tuning post-onset and post-offset for double-peaked neurons

% load('all_peaks.mat')
uCounter = 1;

figure('Position',[1900 500 1800 1200])

colors = colormap(hsv(2));

for pen = 1:length(unit_MIs)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit_MIs{pen,1} '/tmp/Spikes_' unit_MIs{pen,1} '_' unit_MIs{pen,2} '_Good_Pitch.mat']);
    
    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));

%     active = unit_MIs{pen,4};

    for uu = 1:length(units) % for each unit

        unit = units(uu);
        unitSpikes = Y(Y(:,3)==unit,:); % get the spikes of just this unit

        stims_to_plot = stims;

       
        for ss = length(stims):-1:1

            if all_peaks(uCounter,ss,1) == 0 && all_peaks(uCounter,ss,2) == 0
                stims_to_plot(ss) = [];
            end
        end

        if ~isempty(stims_to_plot)

            figure(1); clf
            plot_tuning_by_cond(Y,type,F0,unit,stims_to_plot,zeros(length(stims_to_plot),1),unit_MIs{pen,1},unit_MIs{pen,2},[0 0.1])
    %             hold on
    
            figure(2); clf
            if strcmp(unit_MIs{pen,1},'Noah')
                plot_tuning_by_cond(Y,type,F0,unit,stims_to_plot,zeros(length(stims_to_plot),1),unit_MIs{pen,1},unit_MIs{pen,2},[0.3 0.4])
            else
                plot_tuning_by_cond(Y,type,F0,unit,stims_to_plot,zeros(length(stims_to_plot),1),unit_MIs{pen,1},unit_MIs{pen,2},[0.2 0.3])
            end

%         legend('first','second')

            pause
        end

%         end % ends stimulus loop

        uCounter = uCounter + 1;

    end % ends unit loop
end