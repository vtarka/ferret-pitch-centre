% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, May 2023

figure;
lw = 4;

%% Harmonicity Unit

peak = 12;
width = 2;
f = 1:0.25:17;

tuning_model = exp(-((f-peak)/width).^2/2);

subplot(1,3,1)
hold on
plot(f,tuning_model,'k','linewidth',lw)

peak = 12.2;
width = 1.7;
tuning_model = exp(-((f-peak)/width).^2/2);

plot(f,tuning_model,'b','linewidth',lw)

jitter = rand(17,1)/10;
unr_response = 0.07+jitter;
plot(1:17,unr_response,'r','linewidth',lw)

set(gca,'fontsize',24)

yticks([0 1])
ylabel('Norm. Response')

xlim([1 17])
xticks(1:4:17)
xticklabels(Flist(1:4:17))
xlabel('Pitch (Hz)')

title('Harmonicity Neuron')


%% Temporal Unit

peak = 7;
width = 2.2;
f = 1:0.25:17;

tuning_model = exp(-((f-peak)/width).^2/2);

subplot(1,3,2)
hold on
plot(f,tuning_model,'k','linewidth',lw)

peak = 6.5;
width = 2.4;
tuning_model = exp(-((f-peak)/width).^2/2);

plot(f,tuning_model,'r','linewidth',lw)

set(gca,'fontsize',24)

yticks([0 1])

xlim([1 17])
xticks(1:4:17)
xticklabels(Flist(1:4:17))
xlabel('Pitch (Hz)')

title('Temporal Neuron')

%% Pitch Unit

peak = 10;
width = 2;
f = 1:0.25:17;

tuning_model = exp(-((f-peak)/width).^2/2);

subplot(1,3,3)
hold on
plot(f,tuning_model,'k','linewidth',lw)

peak = 9.9;
width = 2.5;
tuning_model = exp(-((f-peak)/width).^2/2);

plot(f,tuning_model,'b','linewidth',lw)

peak = 10.3;
width = 1.8;
tuning_model = exp(-((f-peak)/width).^2/2);

plot(f,tuning_model,'r','linewidth',lw)

set(gca,'fontsize',24)

yticks([0 1])

xlim([1 17])
xticks(1:4:17)
xticklabels(Flist(1:4:17))
xlabel('Pitch (Hz)')

title('Pitch Neuron')



%% Frequency Response

figure;
hold on

peak = 10;
width = 20;
f = 1:0.25:17;

tuning_model = exp(-((f-peak)/width).^2/2);

plot(f,tuning_model,'k','linewidth',lw)


peak = 15;
width = 2;
f = 1:0.25:17;

tuning_model = exp(-((f-peak)/width).^2/2);

plot(f,tuning_model,'b','linewidth',lw)


peak = 4.5;
width = 6;
f = 1:0.25:17;

tuning_model = exp(-((f-peak)/width).^2/2);

plot(f,tuning_model,'r','linewidth',lw)

set(gca,'fontsize',24)

yticks([0 1])
ylabel('Norm. Response')

xlim([1 17])
xticks(1:4:17)
xticklabels(Flist(1:4:17))
xlabel('Pitch (Hz)')

title('Frequency Tuned Neuron')



%% Response Profile

rp = zeros(5,17);

%%% 1
peak = 13;
width = 2;
f = 1:17;

tuning_model = exp(-((f-peak)/width).^2/2);
rp(1,:) = tuning_model;

%%%% 2
peak = 13.3;
width = 1.2;

tuning_model = exp(-((f-peak)/width).^2/2);
rp(2,:) = tuning_model;

%%%% 3
peak = 12.7;
width = 1.2;

tuning_model = exp(-((f-peak)/width).^2/2);
rp(3,:) = tuning_model;


%%%% 4
peak = 12.9;
width = 2.2;

tuning_model = exp(-((f-peak)/width).^2/2);
rp(4,:) = tuning_model;


%%%% 5
peak = 13.1;
width = 1.7;

tuning_model = exp(-((f-peak)/width).^2/2);
rp(5,:) = tuning_model;


figure;
imagesc(rp)
yticks(1:5)
yticklabels({'High Harm.','Missing F0','Click Train','CT 5% Jitter','CT 10% Jitter'})
xticks(1:4:17)
xticklabels([250 500 1000 2000 4000])
xlabel('F0 (Hz)')
set(gca,'fontsize',24)




