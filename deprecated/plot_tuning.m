%% Plot the tuning curves for each stim type

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13


for ap = 1:length(Animals)

    if ap<9
        continue
    end

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);
    
    stims = {'CT0','CT40','high','low','tone','alt'};%,'rand'};
    stims_legend = {'CT0','CT0','CT40','CT40','high','high','low','low','tone','tone','alt','alt'};
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    units = unique(Y(:,3));
    window = [0 0.15];
    colors = colormap(jet(length(stims)));
    
    figure('Position',[1500 500 1800 1200])
    for uu = 1:length(units)
    
        clf;
        sgtitle(sprintf('Unit # %d, %s, %s',units(uu),Animals{ap},Pens{ap}))
    
        unitSpikes = Y(Y(:,3)==units(uu),:);
        
        % go through each stim we want to plot
        for ss = 1:length(stims)
    
            nSpikes = zeros(length(repeats),length(Flist)); %creates an array of length repList by length Flist
                
            for ff = 1:length(Flist) % do the below for all the frequencies
        
                stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); 
    
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
            spikes = sum(nSpikes,'all');
    
            %%%%%%%%%%%%%%%%%%% plot the PTSH %%%%%%%%%%%%%%%%%%%%%%%%%
            f = polyfit([1:length(Flist)]',mean(nSpikes)',3);
            x1 = linspace(1,length(Flist));
            y1 = polyval(f,x1);
            p = plot(x1,y1); hold on %,'Color',colors(ss,:),'LineWidth',2)
            scatter(1:length(Flist),mean(nSpikes),100,colors(ss,:),'filled')
            set(p,'Color',colors(ss,:),'LineWidth',2)
%             
%             f = fit([1:length(Flist)]',mean(nSpikes)','gauss2');
%             p = plot(f);%,'Color',colors(ss,:),'LineWidth',2)
%             set(p,'Color',colors(ss,:),'LineWidth',2)
%             plot(mean(nSpikes),'Color',colors(ss,:),'LineWidth',2)
% %             errorbar(1:length(Flist),mean(nSpikes),ste(nSpikes)) 
%             hold on 

        
            xticks(1:17)
            xticklabels(Flist)
            ylabel('Evoked Firing Rate (spike/sec)'); % y axis label is the firing rate
            axis tight

    
        end
        legend(stims_legend)
        pause
    end
end