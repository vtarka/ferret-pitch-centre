%% Plot some pure tone stuff
% DEPENDENCIES:
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, May 2023

figure;

file_names = {'HN_units','TN_units','PN_units'};

stims = {'low','high','CT0','tone'};

plot_titles = {'HNs','TNs','PNs'};

window = [0 0.1];

max_sensitivity = 0;
min_sensitivity = 0;

for fi = 1:length(file_names)

    mat_struct = load(file_names{fi});
    mat_cell = struct2cell(mat_struct);
    units_by_rec = mat_cell{1};

    points = zeros(count_units(units_by_rec),4);
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
            points(unit_counter,2) = corr(profile(1,:)',profile(4,:)','rows','complete');
            points(unit_counter,3) = corr(profile(2,:)',profile(3,:)');
            points(unit_counter,4) = corr(profile(2,:)',profile(4,:)');

            subplot(2,3,fi)
            hold on
            scatter(unit_counter,points(unit_counter,1),100,'k','filled','o')
            scatter(unit_counter,points(unit_counter,2),100,'k','filled','square')

            if points(unit_counter,1) > points(unit_counter,2)
                plot([unit_counter unit_counter],points(unit_counter,1:2),'g','linewidth',2)
            else
                plot([unit_counter unit_counter],points(unit_counter,1:2),'r','linewidth',2)
            end

            subplot(2,3,fi+3)
            hold on
            scatter(unit_counter,points(unit_counter,3),100,'k','filled','o')
            scatter(unit_counter,points(unit_counter,4),100,'k','filled','square')

            if points(unit_counter,3) > points(unit_counter,4)
                plot([unit_counter unit_counter],points(unit_counter,3:4),'g','linewidth',2)
            else
                plot([unit_counter unit_counter],points(unit_counter,3:4),'r','linewidth',2)
            end

            unit_counter = unit_counter + 1;

        end
    end % ends recording loop

    subplot(2,3,fi)
    xlim([0 unit_counter])
    ylim([-1 1])
    xticks(10:10:50)
    yticks([-1 0 1])
    set(gca,'fontsize',24)

    subplot(2,3,fi+3)
    xlim([0 unit_counter])
    ylim([-1 1])
    xticks(10:10:50)
    yticks([-1 0 1])
    set(gca,'fontsize',24)

end % ends file loop

subplot(2,3,1)
title({'Harmonicity Neurons','Low Harm.'})

subplot(2,3,2)
title({'Temporal Neurons','Low Harm.'})

subplot(2,3,3)
title({'Pitch Neurons','Low Harm.'})

subplot(2,3,4)
title({'Harmonicity Neurons','High Harm.'})
xlabel('Neuron')
ylabel('Correlation')

subplot(2,3,5)
title({'Temporal Neurons','High Harm.'})

subplot(2,3,6)
title({'Pitch Neurons','High Harm.'})


