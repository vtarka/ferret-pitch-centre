%% Compare the firing rates of harmonicity neurons in response to pure tones vs complex tones
% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, February 2024

file_name = 'HN_units_new_05.mat';

mat_struct = load(file_name);
mat_cell = struct2cell(mat_struct);
units_by_rec = mat_cell{1};

window = [0 0.1];

figure

uCounter = 1;

max_FRs = zeros(count_units(units_by_rec),2);

for pen = 1:length(units_by_rec)


    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/final/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

    units = units_by_rec{pen,3};
%     units = unique(Y(:,3));
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    stims = {'CT0','tone'};

    for uu = 1:length(units)

        unit = units(uu);

        unitSpikes = Y(Y(:,3)==unit,:); % get spikes for just this unit

        maxSpikes = zeros(length(stims),length(Flist));

        % for each stim type
        for ss = 1:length(stims)
                
            meanSpikes = zeros(length(Flist),1);

            % go through each F0 to find BF
            for ff = 1:length(Flist)

                stimNum = find(strcmp(type,stims(ss)) & (F0==Flist(ff))); % unique name for combination of stim type and F0

                if isempty(stimNum) % if this stim type and fo combo wasn't presented
                    continue
                end

                nSpikes = zeros(length(repeats),1); % initialize space to count the number of spikes per trial repeat

                % for each repeat
                for rr = 1:length(repeats)

                    spikeIDXs = unitSpikes(:,4)==stimNum & unitSpikes(:,5)==repeats(rr) & unitSpikes(:,2)>window(1) & unitSpikes(:,2)<window(2);
                    nSpikes(rr) = sum(spikeIDXs);

                end % ends repeat loop

                spikeRates = nSpikes / diff(window);
                maxSpikes(ss,ff) = mean(spikeRates);

%                 meanSpikes(ff) = mean(nSpikes); % average the spike rates across repeats

            end % ends looping through F0s

        end % ends loop through stims

        nexttile;
        hold on

        for i = 1:length(Flist)

            plot(1:2, [maxSpikes(1,i) maxSpikes(2,i)],'k')

        end

        plot(1:2, [mean(maxSpikes(1,:)) mean(maxSpikes(2,:))], 'r','linewidth',3);

        set(gca,'fontsize',16)

        xlim([0 3])
        xticks([])

        [~,I] = max(maxSpikes(1,:));

        max_FRs(uCounter,:) = maxSpikes(:,I);

        uCounter = uCounter + 1;

    end

end

xticks([1 2])

xticklabels(stims)

%%
figure
hold on
for uu = 1:count_units(units_by_rec)

    plot(1:2, max_FRs(uu,:),'k-o')
end

plot(1:2, mean(max_FRs),'r','linewidth',3)

xlim([0 3])
xticks(1:2)
xticklabels(stims)
ylabel('Mean Firing Rate at best F0')

set(gca,'fontsize',16)