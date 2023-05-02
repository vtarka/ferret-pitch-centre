%% Plot click train tuning ordered by pitch salience to test given metric
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, April 2023

mat_struct = load('HNs_99.mat');
mat_cell = struct2cell(mat_struct);
units_by_rec = mat_cell{1};
clear mat_struct mat_cell

stims = {'CT0','CT5','CT10','CT20','CT40'};

nUnits = count_units(units_by_rec);
unit_counter = 1;
CT_profiles = zeros(nUnits,length(stims),17);
pitch_sensitivity = zeros(nUnits,1);

for pen = 1:length(units_by_rec)
    
    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/tmp02/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

    units = units_by_rec{pen,3};
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
%             else
%                 window = [.2 .3];
%             end
%         end

        window = [0 .1];

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

        pitch_sensitivity(unit_counter) = estimate_pitch_salience_sensitivity(profile,0);
        CT_profiles(unit_counter,:,:) = profile;

        unit_counter = unit_counter + 1;
    end
end

[B,I] = sort(pitch_sensitivity);
colors = colormap(hsv(length(stims)));

figure;
for uu = 1:nUnits
    nexttile; hold on;
    for ss = 1:length(stims)
        plot(1:17,squeeze(CT_profiles(I(uu),ss,:)),'Color',colors(ss,:),'linewidth',1.5)
    end
    axis tight
    xticks([])
%     yticks([])
    title(sprintf('%.2f',B(uu)))
end