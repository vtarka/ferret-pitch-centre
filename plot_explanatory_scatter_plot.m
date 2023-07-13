% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, May 2023

colors = [0 0 1; 1 0 0; 0 0 0];
figure;

%% Q2 (-0.33, 0.8)

subplot(2,2,1)

unit = {'Noah','P08',233};

load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit{1} '/tmp02/Spikes_' unit{1} '_' unit{2} '_Good_Pitch.mat']);
plot_tuning_by_cond(Y,type,F0,unit{3},{'low','high','CT0'},[0 0 0],unit{1},unit{2},[0 0.1],colors)
sgtitle('')

xticks(1:4:17)
Flist = unique(F0);
xticklabels(Flist(1:4:17))
xlabel('Pitch (Hz)')

ylim([0 80])
yticks([40 80])
ylabel('Spike Rate')

set(gca,'fontsize',24)


%% Q1 (0.68,0.44)

subplot(2,2,2)

unit = {'Noah','P03',392};

load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit{1} '/tmp02/Spikes_' unit{1} '_' unit{2} '_Good_Pitch.mat']);
plot_tuning_by_cond(Y,type,F0,unit{3},{'low','high','CT0'},[0 0 0],unit{1},unit{2},[0 0.1],colors)
sgtitle('')

xticks(1:4:17)
Flist = unique(F0);
xticklabels(Flist(1:4:17))
xlabel('Pitch (Hz)')

ylim([0 20])
yticks([10 20])
ylabel('Spike Rate')

set(gca,'fontsize',24)


%% Q3 (-0.39,-0.21)

subplot(2,2,3)

unit = {'Noah','P08',203};

load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit{1} '/tmp02/Spikes_' unit{1} '_' unit{2} '_Good_Pitch.mat']);
plot_tuning_by_cond(Y,type,F0,unit{3},{'low','high','CT0'},[0 0 0],unit{1},unit{2},[0 0.1],colors)
sgtitle('')

xticks(1:4:17)
Flist = unique(F0);
xticklabels(Flist(1:4:17))
xlabel('Pitch (Hz)')

ylim([0 20])
yticks([10 20])
ylabel('Spike Rate')

set(gca,'fontsize',24)


%% Q4 (0.67,-0.23)

subplot(2,2,4)

unit = {'Ronnie','P05',448};

load(['/media/veronica/Kat Data/Veronica/pitch_ephys/DansMATLABData/' unit{1} '/tmp02/Spikes_' unit{1} '_' unit{2} '_Good_Pitch.mat']);
plot_tuning_by_cond(Y,type,F0,unit{3},{'low','high','CT0'},[0 0 0],unit{1},unit{2},[0 0.1],colors)
sgtitle('')

xticks(1:4:17)
Flist = unique(F0);
xticklabels(Flist(1:4:17))
xlabel('Pitch (Hz)')

ylim([0 140])
yticks([70 140])
ylabel('Spike Rate')

set(gca,'fontsize',24)


%% Full scatter

figure;
x = [-0.33 0.76 -0.39 0.67];
y = [0.8 0.37 -0.39 -0.23];
scatter(x,y)

scatter(x,y,400,'markerfacecolor',"#77AC30",'linewidth',2.5,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',"#77AC30")

xlim([-1 1])
ylim([-1 1])
xticks([0 1])
yticks([-1 0 1])
grid on
set(gca,'fontsize',24)
