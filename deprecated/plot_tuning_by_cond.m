function a = plot_tuning_by_cond(Y,type,F0,unit,stims,BFs,animal,pen)

    Flist = unique(F0);
    repeats = unique(Y(:,5));
    window = [0 0.15];
    colors = colormap(jet(2));
    
    figure('Position',[1900 500 1800 1200])

    sgtitle(sprintf('%s, %s unit # %d',animal,pen,unit))

    unitSpikes = Y(Y(:,3)==unit,:);
%     plot_handles_for_legend = {};

    % go through each stim we want to plot
    for ss = 1:length(stims)

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

        nSpikes = nSpikes ./ diff(window); % spikes per second
%         spikes = sum(nSpikes,'all');

        %%%%%%%%%%%%%%%%%%% plot the PTSH %%%%%%%%%%%%%%%%%%%%%%%%%
%         errorbar(Flist,mean(nSpikes),ste(nSpikes),'Color',colors(ss,:))
        meanSpikes = mean(nSpikes);
        h = shadedErrorBar(1:17,meanSpikes,ste(nSpikes),{'Color',colors(ss,:)},1);
        h.mainLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
        h.edge(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
        h.edge(2).Annotation.LegendInformation.IconDisplayStyle = 'off';
        hold on
%         if BFs(ss)~= 0
%             p = scatter(BFs(ss),meanSpikes(BFs(ss)),200,'MarkerFaceColor',colors(ss,:));
%             p.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         end
%         set(gca,'XScale','log')
%         f = polyfit([1:length(Flist)]',mean(nSpikes)',3);
%         x1 = linspace(1,length(Flist));
%         y1 = polyval(f,x1);
%         p = plot(x1,y1); hold on %,'Color',colors(ss,:),'LineWidth',2)
%         scatter(1:length(Flist),mean(nSpikes),150,colors(ss,:),'filled')
%         set(p,'Color',colors(ss,:),'LineWidth',2)
%  
%         hold on;
        xticks(1:17)
        xticklabels(num2str(Flist))
        xlabel('F0')
        ylabel('Evoked Firing Rate (spike/sec)'); % y axis label is the firing rate
        set(gca,'fontsize',16)
        axis tight

%         plot_handles_for_legend(end+1) = num2cell(h.mainLine);

    end
   
%     plot_handles_for_legend = cell2mat(plot_handles_for_legend);
%     legend(stims)
    a=1;
end