%% Plot CT0 tuning correlation with low and high tuning, pitch salience is the color map
% DEPENDENCIES:
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, May 2023

figure;

file_names = {'HN_units_f2','TN_units_f2','PN_units_f2'};

stims = {'low','high','tone','CT5','CT10','CT20','CT40'};

plot_titles = {'Harmonicity Neurons','Temporal Neurons','Pitch Neurons'};

window = [0 0.1];

max_sensitivity = 0;
min_sensitivity = 0;

colors = [0 0 1; 1 0 0; 0 0 0];

sensitivity_bins = [-1:0.05:0.95; -0.95:0.05:1]';

for fi = 1:length(file_names) % for each neuronal subtype

    mat_struct = load(file_names{fi});
    mat_cell = struct2cell(mat_struct);
    units_by_rec = mat_cell{1};

    points = zeros(count_units(units_by_rec),3);
    unit_counter = 1;

    for pen = 1:length(units_by_rec)
            
        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/tmp02/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

        units = units_by_rec{pen,3};
        Flist = unique(F0);
        repeats = unique(Y(:,5));

        for uu = 1:length(units)

            unit = units(uu);

            unitSpikes = Y(Y(:,3)==unit,:); % get spikes for just this unit
            profile = zeros(length(stims),17);
           
             % go through each stim we want to plot
            for ss = 1:length(stims)
        
                nSpikes = zeros(length(repeats),length(Flist)); % allocate space to save spiking info for each trial
                    
                for ff = 1:length(Flist) % for each frequency
            
                    stimNum = find(strcmp(type,stims{ss}) & (F0==Flist(ff))); % unique name for combination of stim type and F0
        
                    if isempty(stimNum) % if this stim type and fo combo wasn't presented
                        nSpikes(:,ff) = 0;
                        continue
                    end
            
                    % for each trial of this stim type
                    for rr = 1:length(repeats)
                    
                        spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                        nSpikes(rr,ff) = sum(spikeIDXs);
        
                    end   
                end
        
                nSpikes = nSpikes ./ diff(window); % spikes per second
                meanSpikes = mean(nSpikes); % average across repeats
                profile(ss,:) = meanSpikes;

                if strcmp('low',stims{ss})
                    profile(ss,1:2) = nan;
                end

            end % ends stim loop

            points(unit_counter,1) = corr(profile(1,:)',profile(3,:)','rows','complete');
            points(unit_counter,2) = corr(profile(2,:)',profile(3,:)');

            points(unit_counter,3) = 0; %estimate_pitch_salience_sensitivity(profile(3:end,:));

            unit_counter = unit_counter + 1;

        end
    end % ends recording loop
    
    sensitivity = points(:,3);
    if max(sensitivity) > max_sensitivity
        max_sensitivity = max(sensitivity);
    end

    if min(sensitivity) < min_sensitivity
        min_sensitivity = min(sensitivity);
    end

    [low_corr_sorted, i] = sort(points(:,1));
    sensitivity_sorted = sensitivity(i);

    binned_low_sensitivity = zeros(size(sensitivity_bins,1),1);

    for ii = 1:size(sensitivity_bins,1)
        corr_in_bin = find(low_corr_sorted > sensitivity_bins(ii,1) & low_corr_sorted < sensitivity_bins(ii,2));

        if ~isempty(corr_in_bin)
            sensitivity_in_bin = sensitivity_sorted(corr_in_bin);
            binned_low_sensitivity(ii) = mean(sensitivity_in_bin);
        else
            binned_low_sensitivity(ii) = nan;
        end
    end

    subplot(1,3,fi)

    if fi==1

%         histogram(points(:,1))
        h = histfit(points(:,1));
        X = h(2).XData;
        Y = h(2).YData;
        Y_norm = (Y-min(Y)) / (max(Y)-min(Y));
        area(X,Y_norm,'FaceColor',[1 0 0],'EdgeColor',[1 0 0],'linewidth',2,'facealpha',0.2);
        ax1 = gca;
%         ax1.DataAspectRatio = [1 1 1];
%         ax1.YLim = [0 eps];
        ax1.Color = 'none';
        set(gca,'fontsize',24)
        ax1.XLim = [-1 1];
        ax1.YLim = [0 10];
        xticks([-1 0 1])
        ax1.DataAspectRatio = [1 20 1];

%           pd = fitdist(points(:,1),'Kernel','Width',0.1);
%           x = -1:0.01:1;
%           y = pdf(pd,x);
%           plot(x,y)
%         histogram(points(:,1),-1:0.05:1)
    elseif fi == 2
%         
        h = histfit(points(:,2));
        X = h(2).XData;
        Y = h(2).YData;
        Y_norm = (Y-min(Y)) / (max(Y)-min(Y));
        area(X,Y_norm,'FaceColor',[0 0 1],'EdgeColor',[0 0 1],'linewidth',2,'facealpha',0.2);
        ax1 = gca;
%         ax1.DataAspectRatio = [1 1 1];
%         ax1.XLim = [0 eps];
        ax1.Color = 'none';
        set(gca,'fontsize',24)
        ax1.XLim = [-1 1];
        xticks([-1 0 1])
        ax1.YLim = [0 12];
        ax1.DataAspectRatio = [1 20 1];
%         histogram(points(:,2),-1:0.05:1)

%           pd = fitdist(points(:,2),'Kernel','Width',0.1);
%           x = -1:0.01:1;
%           y = pdf(pd,x);
%           plot(x,y)
    else
%         scatter(points(:,1),points(:,2),200,'markerfacecolor',"#25523B",'linewidth',2.5,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',"#77AC30")
        scatter(points(:,1),points(:,2),200,'markerfacecolor',"#000000",'linewidth',2.5,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',"#000000")
        xlim([-1 1])
        ylim([-1 1])
        xticks([0 1])
        yticks([-1 0 1])
        grid on
    
        set(gca,'fontsize',24)
    end
%     title(plot_titles{fi})
    
%     colormap jet


end % ends file loop



%% Pure tone scatter overlaid

file_names = {'HN_units_f2','TN_units_f2','PN_units_f2'};
% file_names = {'tone_HNs','tone_TNs','tone_PNs'};

stims = {'F0MaskLow','F0MaskHigh','CT0','CT5','CT10','CT20','CT40'};


plot_titles = {'Harmonicity Neurons','Temporal Neurons','Pitch Neurons'};
% plot_titles = {'','',''};

window = [0 0.1];

max_sensitivity = 0;
min_sensitivity = 0;

colors = [0 0 1; 1 0 0; 0 0 0];

sensitivity_bins = [-1:0.05:0.95; -0.95:0.05:1]';

for fi = 1:length(file_names)

    mat_struct = load(file_names{fi});
    mat_cell = struct2cell(mat_struct);
    units_by_rec = mat_cell{1};

    points = zeros(count_units(units_by_rec),3);
    unit_counter = 1;

    for pen = 1:length(units_by_rec)
            
        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/tmp02/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

        units = units_by_rec{pen,3};
        Flist = unique(F0);
        repeats = unique(Y(:,5));

        for uu = 1:length(units)
            unit = units(uu);

            unitSpikes = Y(Y(:,3)==unit,:); % get spikes for just this unit
            profile = zeros(length(stims),17);
           
             % go through each stim we want to plot
            for ss = 1:length(stims)
        
                nSpikes = zeros(length(repeats),length(Flist)); % allocate space to save spiking info for each trial
                    
                for ff = 1:length(Flist) % for each frequency
            
                    stimNum = find(strcmp(type,stims{ss}) & (F0==Flist(ff))); % unique name for combination of stim type and F0
        
                    if isempty(stimNum) % if this stim type and fo combo wasn't presented
                        nSpikes(:,ff) = 0;
                        continue
                    end
            
                    % for each trial of this stim type
                    for rr = 1:length(repeats)
                    
                        spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                        nSpikes(rr,ff) = sum(spikeIDXs);
        
                    end   
                end
        
                nSpikes = nSpikes ./ diff(window); % spikes per second
                meanSpikes = mean(nSpikes); % average across repeats
                profile(ss,:) = meanSpikes;

                if strcmp('F0MaskLow',stims{ss})
                    profile(ss,1:2) = nan;
                end

            end % ends stim loop

            points(unit_counter,1) = corr(profile(1,:)',profile(3,:)','rows','complete');
            points(unit_counter,2) = corr(profile(2,:)',profile(3,:)');

%             points(unit_counter,1)
%             points(unit_counter,2)
%            
%             figure(2)
%             clf
%             plot_tuning_by_cond(Y,type,F0,unit,{'low','high','CT0'},[0 0 0],units_by_rec{pen,1},units_by_rec{pen,2},[0 0.1],colors)

            points(unit_counter,3) = 0; %estimate_pitch_salience_sensitivity(profile(3:end,:));

            unit_counter = unit_counter + 1;

        end
    end % ends recording loop
    
    sensitivity = points(:,3);
    if max(sensitivity) > max_sensitivity
        max_sensitivity = max(sensitivity);
    end

    if min(sensitivity) < min_sensitivity
        min_sensitivity = min(sensitivity);
    end

    [low_corr_sorted, i] = sort(points(:,1));
    sensitivity_sorted = sensitivity(i);

    binned_low_sensitivity = zeros(size(sensitivity_bins,1),1);

    for ii = 1:size(sensitivity_bins,1)
        corr_in_bin = find(low_corr_sorted > sensitivity_bins(ii,1) & low_corr_sorted < sensitivity_bins(ii,2));

        if ~isempty(corr_in_bin)
            sensitivity_in_bin = sensitivity_sorted(corr_in_bin);
            binned_low_sensitivity(ii) = mean(sensitivity_in_bin);
        else
            binned_low_sensitivity(ii) = nan;
        end
    end

    subplot(1,3,fi); hold on
    

%     scatter(points(:,1),points(:,2),100,points(:,3),'filled','linewidth',2.5,'MarkerFaceAlpha',0.2,'MarkerEdgeColor','flat')
%     scatter(points(:,1),points(:,2),200,'markerfacecolor',"#25523B",'linewidth',2.5,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',"#77AC30")
    if fi==1
%         scatter(points(:,1),zeros(length(points(:,1)),1),50,'markerfacecolor',"#25523B",'linewidth',1,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',"#77AC30")
%         scatter(points(:,1),zeros(length(points(:,1)),1),50,'markerfacecolor',"#0000FF",'linewidth',1,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',"#0000FF")
%         histogram(points(:,1))
        figure(4);
        h = histfit(points(:,1));
        X = h(2).XData;
        Y = h(2).YData;
        Y_norm = (Y-min(Y)) / (max(Y)-min(Y));
        figure(1); subplot(1,3,fi); hold on
        area(X,Y_norm,'FaceColor',[1 0 0],'EdgeColor',[1 0 0],'linewidth',2,'facealpha',0.1,'linestyle','--');
        ax1 = gca;

        ax1.Color = 'none';
        set(gca,'fontsize',24)
        ax1.XLim = [-1 1];
        ax1.YLim = [0 10];
        xticks([-1 0 1])
        ax1.DataAspectRatio = [1 20 1];
%         pause
%           subplot(1,3,fi); hold on
%           pd = fitdist(points(:,1),'Kernel','Width',0.1);
%           x = -1:0.01:1;
%           y = pdf(pd,x);
%           plot(x,y)

%         histogram(points(:,1),-1:0.05:1)
    elseif fi == 2
%         scatter(zeros(length(points(:,2)),1),points(:,2),50,'markerfacecolor',"#25523B",'linewidth',1,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',"#77AC30")
%         scatter(zeros(length(points(:,2)),1),points(:,2),50,'markerfacecolor',"#FF0000",'linewidth',1,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',"#FF0000")
        figure(4);
        h = histfit(points(:,2));
        X = h(2).XData;
        Y = h(2).YData;
        Y_norm = (Y-min(Y)) / (max(Y)-min(Y));
        figure(1); subplot(1,3,fi); hold on
        area(X,Y_norm,'FaceColor',[0 0 1],'EdgeColor',[0 0 1],'linewidth',2,'facealpha',0.1,'linestyle','--');
        ax1 = gca;
%         ax1.DataAspectRatio = [1 1 1];
%         ax1.XLim = [0 eps];
        ax1.Color = 'none';
        set(gca,'fontsize',24)
        ax1.XLim = [-1 1];
        xticks([-1 0 1])
        ax1.YLim = [0 12];
        ax1.DataAspectRatio = [1 20 1];
%         histogram(points(:,2),-1:0.05:1)

%           subplot(1,3,fi); hold on
%           pd = fitdist(points(:,2),'Kernel','Width',0.1);
%           x = -1:0.01:1;
%           y = pdf(pd,x);
%           plot(x,y)
    else
%         scatter(points(:,1),points(:,2),200,'markerfacecolor',"#25523B",'linewidth',2.5,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',"#77AC30")
%         figure(2);
        figure(1);
        subplot(1,3,fi); hold on
        scatter(points(:,1),points(:,2),200,'markerfacecolor',"#000000",'linewidth',2.5,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',"#000000",'marker','*')
        xlim([-1 1])
        ylim([-1 1])
        xticks([0 1])
        yticks([-1 0 1])
        grid on
    
        set(gca,'fontsize',24)
    end
    
    set(gca,'fontsize',24)

end % ends file loop



%% Add shaded areas

% shaded_vertices = {[0 -1; 1 -1; 1 1; 0 1], [-1 0; 1 0; 1 1; -1 1],[0 0; 1 0; 1 1; 0 1]};
% 
% for fi = 1:length(file_names)
%     subplot(1,3,fi)
% 
%     vertices = shaded_vertices{fi};
%     patch(vertices(:,1),vertices(:,2),'k','FaceAlpha',0.06,'EdgeColor','none')
% 
% %     clim([-max_sensitivity max_sensitivity])
% % 
% %     if fi == 1
% %         xlabel('CT and low harmonics tuning corr.')
% %         ylabel('CT and high harmonics tuning corr.')
% %     elseif fi == 3
% % %         cb = colorbar;
% % %         cb.Position = cb.Position + [0.02 0 0 0];
% %     end
% end