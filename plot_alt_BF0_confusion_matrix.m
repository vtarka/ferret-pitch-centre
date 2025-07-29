%% Plot confusion matrix of best frequencies for cos phase vs alt phase stim
% Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, February 2024


file_name = '0225TNs';

mat_struct = load(file_name);
mat_cell = struct2cell(mat_struct);
units_by_rec = mat_cell{1};

window = [0 0.1];

confusion_matrix = zeros(17);

extras = [];

all_pairs = [];

ratings = [];

for pen = 1:length(units_by_rec)


    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' units_by_rec{pen,1} '/Batch0324/Spikes_' units_by_rec{pen,1} '_' units_by_rec{pen,2} '_Good_Pitch.mat']);

    units = units_by_rec{pen,3};
%     units = unique(Y(:,3));
    Flist = unique(F0);
    repeats = unique(Y(:,5));
    stims = {'high','alt'};

    high_stim_num = find(strcmp(stims,'high'));
    alt_stim_num = find(strcmp(stims,'alt'));

    for uu = 1:length(units)

        unit = units(uu);

%         if ~sensitivity(uu,high_stim_num) && ~sensitivity(uu,alt_stim_num)
%             continue
%         end

        unitSpikes = Y(Y(:,3)==unit,:); % get spikes for just this unit

        BFs = zeros(1,2);

        

        % for each stim type
        for ss = 1:length(stims)
                
            meanSpikes = zeros(length(Flist),1);
            allSpikes = zeros(length(Flist),length(repeats));

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

                allSpikes(ff,:) = nSpikes;
                meanSpikes(ff) = mean(nSpikes); % average the spike rates across repeats

            end % ends looping through F0s
            
            [max_val,i] = max(meanSpikes); % find the index where the maximum spike rate occurs

            n_instances = length(find(meanSpikes == max_val));

            if n_instances > 1

                i = median(find(meanSpikes == max_val));
%                 i  = i;
            
            end

%             normSpikes = zscore(meanSpikes);
%             max_val = max(normSpikes);
            aboveThresh = find(meanSpikes>max_val*0.75);
            belowThresh = find(meanSpikes<max_val*0.75);
%             minBF0 = min(aboveThresh);
%             maxBF0 = max(aboveThresh);
%             i = median([minBF0 maxBF0]);

%             sigSpikes = meanSpikes;
%             sigSpikes(belowThresh) = 0;
            d = diff(aboveThresh);
            if length(aboveThresh)==1
                title_word = 'good tuning';
                rating = 1;
            elseif length(unique(d))~=1
                title_word = 'bad tuning';
                rating = 0;
            elseif unique(d)~=1
                title_word = 'bad tuning';
                rating = 0;
            else
                title_word = 'good tuning';
                rating = 1;
            end
            com = get_centre_of_mass(allSpikes);
            
            e = std(allSpikes,0,2) / sqrt(size(allSpikes,2));

            clf
            
%             max_val = max(meanSpikes);
            
            hold on;
            xline(com,'r','linewidth',2)
            xline(i,'b','linewidth',2)
            yline(max_val*0.75,'g--')
            shadedErrorBar(1:17,meanSpikes,e,'k')
            legend('COM','BF')
%             plot(com,meanSpikes(floor(com)),'r','marker','*','markersize',30)

            set(gca,'fontsize',14)
            xlabel('F0')
            ylabel('Norm Firing Rate')
            title(title_word)
            axis tight
%             pause

            if length(i)>1
                i = nan;
            end

            if rating == 1 || rating == 0
                BFs(ss) = i;
            else
                BFs(ss) = nan;
            end

        end % ends loop through stims

        all_pairs = [all_pairs; BFs];
        ratings(end+1) = rating;

        if floor(BFs(1))==BFs(1) && floor(BFs(2))==BFs(2)
            confusion_matrix(BFs(1),BFs(2)) = confusion_matrix(BFs(1),BFs(2))+1;
        else
            extras = [extras;BFs(1) BFs(2)];
        end
%         scatter(BFs(1),BFs(2))
    end

end


%% plot results

% figure;
% imagesc(confusion_matrix')
% set(gca,'ydir','normal')
% set(gca,'fontsize',16)
% xlabel('In Phase BF')
% ylabel('Alt. Phase BF')


%% plot as scatter

figure;
hold on

% for i = 1:17
%     for j = 1:17
%         if confusion_matrix(i,j) ~= 0
%             scatter(i,j,confusion_matrix(i,j)*50,'k','filled')
%         end
%     end
% end


set(gca,'fontsize',14)
xticks(1:4:17)
xticklabels(Flist(1:4:17))
xlim([0 18])
ylim([0 18])
yticks(1:4:17)
yticklabels(Flist(1:4:17))
xlabel('In Phase Centre of Mass')
ylabel('Alt Phase Centre of Mass')

plot([0 18],[0,18],':k','linewidth',2)


for i = 1:size(extras,1)
    scatter(extras(i,1),extras(i,2),50,'k','LineWidth',1)
end


function com = get_centre_of_mass(tuning_curve)

idx = 1:size(tuning_curve,1);
tc_product = tuning_curve;

for ff = 1:length(idx)

    tc_product(ff,:) = tuning_curve(ff,:) * idx(ff);

end

total_sum = sum(sum(tuning_curve));
com = sum(sum(tc_product)) / total_sum;

% idx = 1:length(tuning_curve);
% idx_product = tuning_curve' .* idx;
% total_sum = sum(tuning_curve);
% com = sum(idx_product) / total_sum;
end