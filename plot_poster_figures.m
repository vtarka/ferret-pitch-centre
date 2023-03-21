%% Plot specific figures for Oxford Neuroscience Symposium
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, March 2023


figure;

c1 = [0 0 0; 9 220 9; 255 137 0]/255;
c2 = [0 0 0; 0 0 153; 0 204 204; 255 0 255; 255 0 0]/255;

for pen = 1:length(HN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' HN_units{pen,1} '/tmp02/Spikes_' HN_units{pen,1} '_' HN_units{pen,2} '_Good_Pitch.mat']);

    stims1 = {'CT0','low','high'};
    stims2 = {'CT0','CT5','CT10','CT20','CT40'};

%     Flist = unique(F0);
%     repeats = unique(Y(:,5));
%     allUnits = unique(Y(:,3));

    HNUnits = HN_units{pen,3};

    windows = [0 .06; 0.06 .15];

    for uu = 1:size(HNUnits,1)

        if length(find(HNUnits==HNUnits(uu,1)))<2
            if HNUnits(uu,2)==1
    
%               figure
                clf 
                subplot(2,1,1)
                plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims1,zeros(length(stims1),1),HN_units{pen,1},HN_units{pen,2},windows(1,:),c1);
                xlabel('')
                ylabel('')
                xticks([])
                yticks([])
    
                subplot(2,1,2)
                plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims2,zeros(length(stims2),1),HN_units{pen,1},HN_units{pen,2},windows(1,:),c2);
                xlabel('')
                ylabel('')
                xticks([])
                yticks([])
    
                pause
%                 close all
            elseif HNUnits(uu,2)==2
%               figure
                clf 
                subplot(2,1,1)
                plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims1,zeros(length(stims1),1),HN_units{pen,1},HN_units{pen,2},windows(2,:),c1);
                xlabel('')
                ylabel('')
                xticks([])
                yticks([])
    
                subplot(2,1,2)
                plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims2,zeros(length(stims2),1),HN_units{pen,1},HN_units{pen,2},windows(2,:),c2);
                xlabel('')
                ylabel('')
                xticks([])
                yticks([])
    
                pause
%                 close all
            end
        else

%             figure
            clf
            subplot(2,2,1)
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims1,zeros(length(stims1),1),HN_units{pen,1},HN_units{pen,2},windows(1,:),c1);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

            subplot(2,2,2)
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims2,zeros(length(stims2),1),HN_units{pen,1},HN_units{pen,2},windows(1,:),c2);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

            subplot(2,2,3)
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims1,zeros(length(stims1),1),HN_units{pen,1},HN_units{pen,2},windows(2,:),c1);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

            subplot(2,2,4)
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims2,zeros(length(stims2),1),HN_units{pen,1},HN_units{pen,2},windows(2,:),c2);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

            pause
%             close all
        end

    end
end % ends recording loop


%% Plot HN locations

loc_frequencies = zeros(5,1);
loc_total_units = zeros(5,1);

locs = [1 1 2 1 5 3 5 4 3 3 1 2 3 2 1 1 1 3 4 3];

n = 0;
nU = 0;

for pen = 1:length(HN_units)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' HN_units{pen,1} '/tmp02/Spikes_' HN_units{pen,1} '_' HN_units{pen,2} '_Good_Pitch.mat']);
    
    HNUs = HN_units{pen,3};

    if ~isempty(HNUs)
        uHNUs = unique(HNUs);

        loc_frequencies(locs(pen)) = loc_frequencies(locs(pen)) + length(uHNUs);

    end

    loc_total_units(locs(pen)) = loc_total_units(locs(pen)) + length(unique(Y(:,3)));

end

figure;
H = bar(loc_frequencies./loc_total_units,'k');
xlim([0.5 4.5])
xticks(1:4)
xticklabels({'low A1','high A1','low AAF','high AAF'})
yticks(0:0.05:.25)
yticklabels(0:0.05*100:.25*100)
ylabel('% Harmonicity Neurons')
set(gca,'fontsize',24)


%%
figure;

c3 = [0 246 148; 34 0 255; 246 164 0]/255;

for pen = 1:length(HN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' HN_units{pen,1} '/tmp02/Spikes_' HN_units{pen,1} '_' HN_units{pen,2} '_Good_Pitch.mat']);

    stims1 = {'CT0','high'};
    stims2 = {'CT0','CT5','CT10','CT20','CT40'};
    stims3 = {'alt','rand','high'};

    HNUnits = HN_units{pen,3};

    windows = [0 .06; 0.06 .15];

    for uu = 1:size(HNUnits,1)

        if length(find(HNUnits==HNUnits(uu,1)))<2
            if HNUnits(uu,2)==1
    
%                 figure
                clf 
                subplot(3,1,1)
                plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims1,zeros(length(stims1),1),HN_units{pen,1},HN_units{pen,2},windows(1,:),c1);
                xlabel('')
                ylabel('')
                xticks([])
                yticks([])
    
                subplot(3,1,2)
                plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims2,zeros(length(stims2),1),HN_units{pen,1},HN_units{pen,2},windows(1,:),c2);
                xlabel('')
                ylabel('')
                xticks([])
                yticks([])

                subplot(3,1,3)
                plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims3,zeros(length(stims3),1),HN_units{pen,1},HN_units{pen,2},windows(1,:),c3);
                xlabel('')
                ylabel('')
                xticks([])
                yticks([])
    
                pause
%                 close all
            elseif HNUnits(uu,2)==2
%                 figure
                clf 
                subplot(3,1,1)
                plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims1,zeros(length(stims1),1),HN_units{pen,1},HN_units{pen,2},windows(2,:),c1);
                xlabel('')
                ylabel('')
                xticks([])
                yticks([])
    
                subplot(3,1,2)
                plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims2,zeros(length(stims2),1),HN_units{pen,1},HN_units{pen,2},windows(2,:),c2);
                xlabel('')
                ylabel('')
                xticks([])
                yticks([])

                subplot(3,1,3)
                plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims3,zeros(length(stims3),1),HN_units{pen,1},HN_units{pen,2},windows(2,:),c3);
                xlabel('')
                ylabel('')
                xticks([])
                yticks([])

                pause
%                 close all
            end
        else

%             figure
            clf 
            subplot(3,2,1)
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims1,zeros(length(stims1),1),HN_units{pen,1},HN_units{pen,2},windows(1,:),c1);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

            subplot(3,2,3)
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims2,zeros(length(stims2),1),HN_units{pen,1},HN_units{pen,2},windows(1,:),c2);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])


            subplot(3,2,5)
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims3,zeros(length(stims3),1),HN_units{pen,1},HN_units{pen,2},windows(1,:),c3);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

            subplot(3,2,2)
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims1,zeros(length(stims1),1),HN_units{pen,1},HN_units{pen,2},windows(2,:),c1);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

            subplot(3,2,4)
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims2,zeros(length(stims2),1),HN_units{pen,1},HN_units{pen,2},windows(2,:),c2);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])


            subplot(3,2,6)
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims3,zeros(length(stims3),1),HN_units{pen,1},HN_units{pen,2},windows(2,:),c3);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])


            pause
%             close all
        end

    end
end % ends recording loop


%%
figure;

c = colormap(hsv(3));

for pen = 1:length(HN_units)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' HN_units{pen,1} '/tmp02/Spikes_' HN_units{pen,1} '_' HN_units{pen,2} '_Good_Pitch.mat']);

    stims1 = {'CT0','low','high'};
    stims2 = {'CT0','CT5','CT10','CT20','CT40'};
%     stims3 = {'alt','rand','high'};

    HNUnits = HN_units{pen,3};

    windows = [0 .06; 0.06 .15];

    for uu = 1:size(HNUnits,1)

        if length(find(HNUnits==HNUnits(uu,1)))<2
            if HNUnits(uu,2)==1
    
%                 figure
                clf 
                subplot(3,1,1)
                plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims1,zeros(length(stims1),1),HN_units{pen,1},HN_units{pen,2},windows(1,:),c);
                xlabel('')
                ylabel('')
                xticks([])
                yticks([])
    
                subplot(3,1,2)
                plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims2,zeros(length(stims2),1),HN_units{pen,1},HN_units{pen,2},windows(1,:));
                xlabel('')
                ylabel('')
                xticks([])
                yticks([])
    
                pause
%                 close all
            elseif HNUnits(uu,2)==2
%                 figure
                clf 
                subplot(3,1,1)
                plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims1,zeros(length(stims1),1),HN_units{pen,1},HN_units{pen,2},windows(2,:),c);
                xlabel('')
                ylabel('')
                xticks([])
                yticks([])
    
                subplot(3,1,2)
                plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims2,zeros(length(stims2),1),HN_units{pen,1},HN_units{pen,2},windows(2,:));
                xlabel('')
                ylabel('')
                xticks([])
                yticks([])

                pause
%                 close all
            end
        else

%             figure
            clf 
            subplot(3,2,1)
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims1,zeros(length(stims1),1),HN_units{pen,1},HN_units{pen,2},windows(1,:),c);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

            subplot(3,2,3)
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims2,zeros(length(stims2),1),HN_units{pen,1},HN_units{pen,2},windows(1,:));
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

            subplot(3,2,2)
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims1,zeros(length(stims1),1),HN_units{pen,1},HN_units{pen,2},windows(2,:),c);
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])

            subplot(3,2,4)
            plot_tuning_by_cond(Y,type,F0,HNUnits(uu,1),stims2,zeros(length(stims2),1),HN_units{pen,1},HN_units{pen,2},windows(2,:));
            xlabel('')
            ylabel('')
            xticks([])
            yticks([])



            pause
%             close all
        end

    end
end % ends recording loop