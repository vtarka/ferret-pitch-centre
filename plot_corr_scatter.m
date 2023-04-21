%% Plot a scatter plot of tuning correlation vs pitch sensitivity for the 3 populations of neurons
% DEPENDENCIES:
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, April 2023

files = {'HNs_fNoah','TNs_fNoah','PNs_fNoah'};
plot_titles = {'Harmonicity Neurons','Temporal Neurons','Pitch Neurons'};

binary_yn = 'n'; % set this to y if you want the pitch sensitivity to be shown as either yes or no, rather than a continuous metric

stims = {'CT0','CT5','CT10','CT20','CT40','high','low'};

for ff = 1:length(files)

    units_by_rec = load(files{ff});

    nUnits = count_units(units_by_rec);
    corr_v_sens = zeros(nUnits,4);
    unit_counter = 1;

    for pen = 1:length(units_by_rec)
        
        load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/tmp02/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

        units = units_by_rec{pen,3};

        if ~isempty(units)
            units(units(:,2)==3,:) = [];
            units = unique(units(:,1));
        else
            continue
        end

        Flist = unique(F0);
        repeats = unique(Y(:,5));

        for uu = 1:size(units,1)

            unit = units(uu,1);

            unitSpikes = Y(Y(:,3)==unit,:); % get spikes for just this unit
            profile = zeros(length(stims),17);

    %         if units(uu,2)==1
    %             window = [0 0.06];
    %         elseif units(uu,2)==2
    %             window = [0.06 0.15];
    %         else
    %             if pen<9
    %                 window = [.3 .4];
    %             elsef
    %                 window = [.2 .3];
    %             end
    %         end

            window = [0 0.085];

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

    %         [r_low,~,~] = bootstrap_corr(Y,type,F0,{'CT0','low'},unit,window,100);
    %         [r_high,~,~] = bootstrap_corr(Y,type,F0,{'CT0','high'},unit,window,100);
    %         corr_v_sens(unit_counter,1) = r_low;
    %         corr_v_sens(unit_counter,2) = r_high;

            corr_v_sens(unit_counter,1) = corr(profile(1,:)',profile(7,:)','rows','complete');
            corr_v_sens(unit_counter,2) = corr(profile(1,:)',profile(6,:)');
            corr_v_sens(unit_counter,3) = estimate_pitch_salience_sensitivity(profile(1:5,:));
            corr_v_sens(unit_counter,4) = corr(profile(1,:)',max(profile(6:7,:))');

            unit_counter = unit_counter + 1;

        end
    end

    if strcmp(binary_yn,'y')
        r = [255 0 0];
        b = [0 0 0];
        colors = zeros(length(corr_v_sens),3);
        colors(corr_v_sens(:,3)==1,1) = r(1);
    else
        colors = corr_v_sens(:,3);
    end

    figure(2); subplot(1,3,ff)
    scatter(corr_v_sens(:,1),corr_v_sens(:,2),150, colors,'filled',"o")
    xlim([-1 1])
    ylim([-1 1])
    xticks([-1 0 1])
    yticks([-1 0 1])

    title(plot_titles{ff})

    set(gca,'fontsize',24)
end


%%

% x = corr_v_sens(:,1); y = corr_v_sens(:,3);
% p = polyfit(x,y,1);
% xFit = linspace(min(x),max(x),1000);
% yFit = polyval(p,xFit);
% plot(x,y,'k.','markersize',25)
% plot(xFit, yFit,'k-','linewidth',2)
% 
% x = corr_v_sens(:,2);
% p = polyfit(x,y,1);
% xFit = linspace(min(x),max(x),1000);
% yFit = polyval(p,xFit);
% plot(x,y,'r.','markersize',25)
% plot(xFit, yFit,'r-','linewidth',2)

% x = corr_v_sens(:,4);
% p = polyfit(x,y,1);
% xFit = linspace(min(x),max(x),1000);
% yFit = polyval(p,xFit);
% plot(x,y,'b.','markersize',25)
% plot(xFit, yFit,'b-','linewidth',2)