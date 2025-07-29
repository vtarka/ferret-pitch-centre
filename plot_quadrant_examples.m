%% Plot example neurons from each quadrant of the summary scatter plots
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, March 2024

load('CT_points.mat')
CT_points = points;
load('tone_points.mat')
tone_points = points;
clear points

stims = {'low','high','CT0','tone'};

load('PN_units')

window = [0 0.1];

colors = [0 0 1; 1 0 0; 0 0 0; 0 0 0];

HL_colors = [0.4660 0.6740 0.1880; 0.4940 0.1840 0.5560; 0.8500 0.3250 0.0980];
color_titles = {'green','purple','orange'};

uCounter = 1;
fig_counter = 1;
figure;

for pen = 1:length(PN_units)
            
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' PN_units{pen,1} '/Batch0324/Spikes_' PN_units{pen,1} '_' PN_units{pen,2} '_Good_Pitch.mat']);

    units = PN_units{pen,3};
    Flist = unique(F0);
    repeats = unique(Y(:,5));

    for uu = 1:length(units)
        unit = units(uu);

        if ismember(uCounter,[2 3 17])
            figure(3);
            subplot(1,3,fig_counter)
%             subplot(1,2,1)
            hold off
            plot_tuning_by_cond(Y,type,F0,unit,stims,[0 0 0 0],PN_units{pen,1},PN_units{pen,2},window,colors);
            title(color_titles{fig_counter})
            set(gca,'fontsize',20)
            xticks(1:4:17)
            xticklabels(Flist(1:4:17))
            
            
            figure (1)
            subplot(1,3,3)
%             subplot(1,2,2)
%             scatter([CT_points(uCounter,1) tone_points(uCounter,1)],[CT_points(uCounter,2) tone_points(uCounter,2)],50,[1 0 1; 0 1 0],'filled')
            scatter(tone_points(uCounter,1),tone_points(uCounter,2),200,'markerfacecolor',HL_colors(fig_counter,:),'linewidth',2.5,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',HL_colors(fig_counter,:))
%             ylim([-1 1])
%             xlim([-1 1])
%             xline(0)
%             yline(0)
    
            figure (2)
            subplot(1,3,3)
            scatter(CT_points(uCounter,1),CT_points(uCounter,2),200,'markerfacecolor',HL_colors(fig_counter,:),'linewidth',2.5,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',HL_colors(fig_counter,:))

           fig_counter = fig_counter + 1;

        end
        uCounter = uCounter + 1;
        


    end

% 
%         for uu = 1:length(units)
%         unit = units(uu);
% 
% %         if ismember(uCounter,[3 9 17])
%             figure(3);
% %             subplot(1,3,fig_counter)
%             subplot(1,2,1)
%             hold off
%             plot_tuning_by_cond(Y,type,F0,unit,stims,[0 0 0 0],PN_units{pen,1},PN_units{pen,2},window,colors);
% %             title(color_titles{fig_counter})
%             
%             
%             
% %             figure (1)
% %             subplot(1,3,3)
%             subplot(1,2,2)
%             scatter([CT_points(uCounter,1) tone_points(uCounter,1)],[CT_points(uCounter,2) tone_points(uCounter,2)],50,[1 0 1; 0 1 0],'filled')
% %             scatter(tone_points(uCounter,1),tone_points(uCounter,2),200,'markerfacecolor',HL_colors(fig_counter,:),'linewidth',2.5,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',HL_colors(fig_counter,:))
%             ylim([-1 1])
%             xlim([-1 1])
%             xline(0)
%             yline(0)
%     
% %             figure (2)
% %             subplot(1,3,3)
% %             scatter(CT_points(uCounter,1),CT_points(uCounter,2),200,'markerfacecolor',HL_colors(fig_counter,:),'linewidth',2.5,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',HL_colors(fig_counter,:))
% % 
% %            fig_counter = fig_counter + 1;
% 
% %         end
%         uCounter = uCounter + 1;
% 
%         pause
%         
% 
% 
%     end
end