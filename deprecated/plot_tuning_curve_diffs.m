%% Plot
% For each unit, take the absolute difference between every possible
% pairing of stim types

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

figure('Position',[2100 500 1600 1000])



for ap = 1:length(Animals)

    if ap<9
        continue
    end

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    
    stims = unique(type);
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));
    window = [0 0.15];

    % for each unit
    for uu = 1:length(units)
        unit = units(uu); % get this unit's number

        % if this unit wasn't F0-sensitive for any condition, skip it
        if isempty(find(sensitivity(uu,:), 1))
            continue
        end

        % get our figure set up
        clf;
        sgtitle(sprintf('%s, %s unit # %d',Animals{ap},Pens{ap},unit))

        unitSpikes = Y(Y(:,3)==unit,:); % get the spikes of just this unit

        diffMatrix = zeros(length(stims));

        for ss = 1:length(stims)

            nSpikes = zeros(length(repeats),length(Flist));

            for ff = 1:length(Flist)
    
                stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); 
    
                if isempty(stimNum) % if this stim type and fo combo wasn't presented
                    continue
                end
            
                for rr = 1:length(repeats)

                    spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                    nSpikes(rr,ff) = sum(spikeIDXs);

                end % ends repeat loop
            end % ends F0 loop

            meanSpikes_s1 = mean(nSpikes);
            
            s=1;
            while s < ss
%             for s = 1:length(stims)

%                 if ss==s
%                     diffMatrix(ss,s) = 0;
%                     continue
%                 end
              

                nSpikes = zeros(length(repeats),length(Flist));
    
                for ff = 1:length(Flist)
        
                    stimNum = find(strcmp(type,stims(s)) & (F0==Flist(ff))); 
        
                    if isempty(stimNum) % if this stim type and fo combo wasn't presented
                        continue
                    end
                
                    for rr = 1:length(repeats)
    
                        spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                        nSpikes(rr,ff) = sum(spikeIDXs);
    
                    end % ends repeat loop
                end % ends F0 loop
    
                meanSpikes_s2 = mean(nSpikes);

                diffMatrix(ss,s) = sum(abs(diff(meanSpikes_s1 - meanSpikes_s2)));
                s = s+1;

            end % ends inner stim loop

        end % ends stim loop

        diffMatrix = diffMatrix;
%         p = symrcm(diffMatrix);
%         reorderedMatrix = diffMatrix(p,p);
        imagesc(diffMatrix); colorbar; colormap(flipud(parula));
        xticks(1:13); xticklabels(stims(p))
        yticks(1:13); yticklabels(stims(p))
        pause
%         cb = colorbar;
%         cb.Ticks = floor(min(min(diffMatrix))):ceil(max(max(diffMatrix)));
%         set(gca,'YDir','reverse')
%         xticks(1:17); xticklabels(Flist)
%         yticks(1:length(stims_to_plot)); yticklabels(stims_to_plot)
%         hold on; 
%         colormap default;
%         
%         for bf = 1:length(BFs_to_plot)
%             plot(BFs_to_plot(bf),bf,'r+', 'LineWidth',10,'MarkerSize',30)
%         end % ends loop through BFs 
% 
%         plot_tuning_by_cond(Y,type,F0,unit,stims_to_plot,BFs(uu,stim_to_plot_IDX),Animals{ap},Pens{ap});
%         figure(1); colormap default;
%         pause
%         close 2

    end % ends unit loop
end % ends recording loop