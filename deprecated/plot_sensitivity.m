
%% Plot how many neurons are F0-sensitive to each stimulus

Animals = {'Noah','Noah','Noah','Noah','Noah','Noah','Noah','Noah',...
    'Ronnie','Ronnie','Ronnie','Ronnie','Derry','Derry','Derry','Derry',...
    'Dory','Dory','Dory','Dory'};

Pens = {'P01','P02','P03','P04','P05','P06','P07','P08',...
    'P04','P05','P08','P13','P02','P03','P05','P08',...
    'P00','P01','P02','P04'};

Qualia = 'Good';

% %stimList: 'CT0'    'CT10'    'CT20'    'CT40'    'CT5'    'F0MaskHigh'    'F0MaskLow'    'allHarm'      'alt'     'high'    'low'    'rand'    'tone'
% %             1       2          3         4        5             6          7                 8           9          10       11       12        13

nSensitive = zeros(length(Animals),13);

for ap = 1:length(Animals)

    load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' Animals{ap} '/tmp/Spikes_' Animals{ap} '_' Pens{ap} '_Good_Pitch.mat']);

    stims = unique(type);
%     Flist = unique(F0);

    for ss = 1:length(stims)

        nSensitive(ap,ss) = length(find(BFs(:,ss)));

    end

end


% plot all the penetrations grouped together
figure; histogram('Categories',string(stims),'BinCounts',sum(nSensitive))

% separate the data by penetration
h = [];
for pp = 1:length(Pens)
     h = [h; nSensitive(pp,:)];
end

figure;
colors = colormap(jet(20));
b = bar(1:13,h','FaceColor','flat');
for bb = 1:size(h',2)
    b(bb).CData = colors(bb,:);
end
legend(Pens)
xticklabels(stims)
set(gca,'FontSize',18)

% flip the previous
h = [];
for pp = 1:length(Pens)
     h = [h; nSensitive(pp,:)];
end

figure;
subplot(2,1,1)
colors = colormap(jet(13));
b = bar(1:10,h(1:10,:),'FaceColor','flat');
for bb = 1:size(h,2)
    b(bb).CData = colors(bb,:);
end
xticklabels({'N1','N2','N3','N4','N5','N6','N7','N8','R4','R5'});
set(gca,'FontSize',16)
xlabel('Penetration')
ylabel('# of Units F0-sensitive')
legend(stims,'Location','northwest')

subplot(2,1,2)
colors = colormap(jet(13));
b = bar(1:10,h(11:20,:),'FaceColor','flat');
for bb = 1:size(h,2)
    b(bb).CData = colors(bb,:);
end
xticklabels({'R8','R13','De2','De3','De5','De8','Do0','Do1','Do2','Do4'})
set(gca,'FontSize',16)
xlabel('Penetration')
ylabel('# of Units F0-sensitive')

% separate the data by animal
h = [sum(nSensitive(1:8,:)); sum(nSensitive(9:12,:)); sum(nSensitive(13:16,:)); sum(nSensitive(17:20,:))];

figure;
b = bar(1:13,h','FaceColor','flat');
colors = colormap(jet(4));
for bb = 1:size(h',2)
    b(bb).CData = colors(bb,:);
end
legend('Noah','Ronnie','Derry','Dory')
xticklabels(stims)


% group each stim by animal
figure;
b = bar(1:4,h,'FaceColor','flat');
colors = colormap(jet(13));
for bb = 1:size(h,2)
    b(bb).CData = colors(bb,:);
end
legend(stims)
xticklabels({'Noah','Ronnie','Derry','Dory'})
set(gca,'FontSize',20)
xlabel('Ferret')
ylabel('# of Units F0-sensitive')


% group each stim by area
areas = {'low A1','high A1','high AAF','low AAF','low PPF'};
h = [sum(nSensitive([1,2,4,11,15,16,17],:)); sum(nSensitive([3,14,12],:));...
    sum(nSensitive([8,19],:)); sum(nSensitive([9,10,6,13,20,18],:));...
    sum(nSensitive([5,7],:))];

figure;
b = bar(1:5,h,'FaceColor','flat');
colors = colormap(jet(13));
for bb = 1:size(h,2)
    b(bb).CData = colors(bb,:);
end
legend(stims)
xticklabels(areas)
set(gca,'FontSize',20)
ylabel('# of Units F0-sensitive')



%% Normalizing everything

% plot all the penetrations grouped together
figure; histogram('Categories',string(stims),'BinCounts',sum(nSensitive))

penUnitCounts = [27,146,156,211,36,46,9,43,...
    34,77,146,58,60,61,37,29,...
    43,40,26,26];

% separate the data by penetration
h = [];
for pp = 1:length(Pens)
     h = [h; nSensitive(pp,:)];
end

figure;
colors = colormap(jet(20));
b = bar(1:13,h','FaceColor','flat');
for bb = 1:size(h',2)
    b(bb).CData = colors(bb,:);
end
legend(Pens)
xticklabels(stims)
set(gca,'FontSize',18)

% flip the previous
h = [];
for pp = 1:length(Pens)
     h = [h; nSensitive(pp,:)./penUnitCounts(pp).*100];
end

figure;
subplot(2,1,1)
colors = colormap(jet(13));
b = bar(1:10,h(1:10,:),'FaceColor','flat');
for bb = 1:size(h,2)
    b(bb).CData = colors(bb,:);
end
xticklabels({'N1','N2','N3','N4','N5','N6','N7','N8','R4','R5'});
set(gca,'FontSize',16)
xlabel('Penetration')
ylabel('% of Units F0-sensitive')
ylim([0 90])
legend(stims,'Location','northwest')

subplot(2,1,2)
colors = colormap(jet(13));
b = bar(1:10,h(11:20,:),'FaceColor','flat');
for bb = 1:size(h,2)
    b(bb).CData = colors(bb,:);
end
xticklabels({'R8','R13','De2','De3','De5','De8','Do0','Do1','Do2','Do4'})
set(gca,'FontSize',16)
xlabel('Penetration')
ylabel('% of Units F0-sensitive')
ylim([0 90])

% separate the data by animal
h = [sum(nSensitive(1:8,:))./sum(penUnitCounts(1:8))*100;...
    sum(nSensitive(9:12,:))./sum(penUnitCounts(9:12))*100;...
    sum(nSensitive(13:16,:))./sum(penUnitCounts(13:16))*100;...
    sum(nSensitive(17:20,:))./sum(penUnitCounts(17:20))*100];

figure;
b = bar(1:4,h,'FaceColor','flat');
colors = colormap(jet(13));
for bb = 1:size(h,2)
    b(bb).CData = colors(bb,:);
end
legend(stims)
xticklabels({'Noah','Ronnie','Derry','Dory'})
set(gca,'FontSize',20)
xlabel('Ferret')
ylabel('% of Units F0-sensitive')


% group each stim by area
areas = {'low A1','high A1','high AAF','low AAF','low PPF'};
h = [sum(nSensitive([1,2,4,11,15,16,17],:))./sum(penUnitCounts([1,2,4,11,15,16,17]))*100;...
    sum(nSensitive([3,14,12],:))./sum(penUnitCounts([3,14,12]))*100;...
    sum(nSensitive([8,19],:))./sum(penUnitCounts([8,19]))*100;...
    sum(nSensitive([9,10,6,13,20,18],:))./sum(penUnitCounts([9,10,6,13,20,18]))*100;...
    sum(nSensitive([5,7],:))./sum(penUnitCounts([5,7]))*100];

figure;
b = bar(1:5,h,'FaceColor','flat');
colors = colormap(jet(13));
for bb = 1:size(h,2)
    b(bb).CData = colors(bb,:);
end
legend(stims)
xticklabels(areas)
set(gca,'FontSize',20)
ylabel('# of Units F0-sensitive')
