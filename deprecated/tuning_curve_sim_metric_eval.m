%% Veronica Tarka
% January 2023
% veronica.tarka@dpag.ox.ac.uk

%% Testing out toy examples of tuning curves to see which metrics assess their similarity the most accurately


%% Toy example 1: identical tuning

figure; 

R = zeros(4,1);
BF_diff = zeros(4,1);
ABS_diff = zeros(4,1);
BW_diff = zeros(4,1);

t1 = [0,0,0,1,2,1,0,0,0,0,0,0,0,0,0,0,0]';
t2 = t1;

corr_estmt = corr(t1,t2);
R(1) = corr_estmt;

[~,bf1] = max(t1);
[~,bf2] = max(t2);
bf_diff_estmt = abs(bf1 - bf2);
BF_diff(1) = bf_diff_estmt;

abs_diff_estmt = sum(abs(t1-t2));
ABS_diff(1) = abs_diff_estmt;

t1_gauss = fit([1:17]',t1,'gauss1');
t1_bw = 2*t1_gauss.c1*sqrt(log(2));

t2_gauss = fit([1:17]',t2,'gauss1');
t2_bw = 2*t2_gauss.c1*sqrt(log(2));
bw_diff_estmt = abs(t1_bw - t2_bw);
BW_diff(1) = bw_diff_estmt;

subplot(4,2,1)
plot(1:17,t1,'k','LineWidth',3); hold on;
plot(1:17,t2,'--r','LineWidth',2); 
xlim([0 19])
% text(10,max(t1)-1,sprintf('R: %.2f\nBF diff: %d\nABS diff: %d\nBW diff: %0.2f',corr_estmt,bf_diff_estmt,abs_diff_estmt,bw_diff_estmt),'FontSize',18)


% slightly shifted

t2 = [0,0,0,0,0,1,2,1,0,0,0,0,0,0,0,0,0]';

corr_estmt = corr(t1,t2);
R(2) = corr_estmt;

[~,bf1] = max(t1);
[~,bf2] = max(t2);
bf_diff_estmt = abs(bf1 - bf2);
BF_diff(2) = bf_diff_estmt;

abs_diff_estmt = sum(abs(t1-t2));
ABS_diff(2) = abs_diff_estmt;

t1_gauss = fit([1:17]',t1,'gauss1');
t1_bw = 2*t1_gauss.c1*sqrt(log(2));

t2_gauss = fit([1:17]',t2,'gauss1');
t2_bw = 2*t2_gauss.c1*sqrt(log(2));
bw_diff_estmt = abs(t1_bw - t2_bw);
BW_diff(2) = bw_diff_estmt;

subplot(4,2,3)
plot(1:17,t1,'k','LineWidth',3); hold on;
plot(1:17,t2,'--r','LineWidth',2); 
xlim([0 19])
% text(10,max(t1)-1,sprintf('R: %.2f\nBF diff: %d\nABS diff: %d\nBW diff: %0.2f',corr_estmt,bf_diff_estmt,abs_diff_estmt,bw_diff_estmt),'FontSize',18)



% more shifted

t2 = [0,0,0,0,0,0,0,0,1,2,1,0,0,0,0,0,0]';

corr_estmt = corr(t1,t2);
R(3) = corr_estmt;

[~,bf1] = max(t1);
[~,bf2] = max(t2);
bf_diff_estmt = abs(bf1 - bf2);
BF_diff(3) = bf_diff_estmt;

abs_diff_estmt = sum(abs(t1-t2));
ABS_diff(3) = abs_diff_estmt;

t1_gauss = fit([1:17]',t1,'gauss1');
t1_bw = 2*t1_gauss.c1*sqrt(log(2));

t2_gauss = fit([1:17]',t2,'gauss1');
t2_bw = 2*t2_gauss.c1*sqrt(log(2));
bw_diff_estmt = abs(t1_bw - t2_bw);
BW_diff(3) = bw_diff_estmt;

subplot(4,2,5)
plot(1:17,t1,'k','LineWidth',3); hold on;
plot(1:17,t2,'--r','LineWidth',2); 
xlim([0 19])
% text(10,max(t1)-1,sprintf('R: %.2f\nBF diff: %d\nABS diff: %d\nBW diff: %0.2f',corr_estmt,bf_diff_estmt,abs_diff_estmt,bw_diff_estmt),'FontSize',18)


% most shifted

t2 = [0,0,0,0,0,0,0,0,0,0,0,0,1,2,1,0,0]';

corr_estmt = corr(t1,t2);
R(4) = corr_estmt;

[~,bf1] = max(t1);
[~,bf2] = max(t2);
bf_diff_estmt = abs(bf1 - bf2);
BF_diff(4) = bf_diff_estmt;

abs_diff_estmt = sum(abs(t1-t2));
ABS_diff(4) = abs_diff_estmt;

t1_gauss = fit([1:17]',t1,'gauss1');
t1_bw = 2*t1_gauss.c1*sqrt(log(2));

t2_gauss = fit([1:17]',t2,'gauss1');
t2_bw = 2*t2_gauss.c1*sqrt(log(2));
bw_diff_estmt = abs(t1_bw - t2_bw);
BW_diff(4) = bw_diff_estmt;

subplot(4,2,7)
plot(1:17,t1,'k','LineWidth',3); hold on;
plot(1:17,t2,'--r','LineWidth',2); 
xlim([0 19])
% text(10,max(t1)-1,sprintf('R: %.2f\nBF diff: %d\nABS diff: %d\nBW diff: %0.2f',corr_estmt,bf_diff_estmt,abs_diff_estmt,bw_diff_estmt),'FontSize',18)



% plot the metrics for each of these curves

subplot(4,2,[2 4 6 8]); hold on;
plot(1:4,1-R,'LineWidth',2)
plot(1:4,BF_diff,'LineWidth',2)
plot(1:4,ABS_diff,'LineWidth',2)
plot(1:4, BW_diff,'LineWidth',2)

legend({'1-Corr','BF diff','Absolute diff','Bandwidth diff'},'Location','northwest')
set(gca,'Fontsize',16)

% set(gca,'YDir','reverse')

%% Toy example 2: differing bandwidth

figure; 

R = zeros(4,1);
BF_diff = zeros(4,1);
ABS_diff = zeros(4,1);
BW_diff = zeros(4,1);

t1 = [0,0,0,1,2,1,0,0,0,0,0,0,0,0,0,0,0]';
t2 = [0,0.5,1,1,2,1,1,0.5,0,0,0,0,0,0,0,0,0]';

corr_estmt = corr(t1,t2);
R(1) = corr_estmt;

[~,bf1] = max(t1);
[~,bf2] = max(t2);
bf_diff_estmt = abs(bf1 - bf2);
BF_diff(1) = bf_diff_estmt;

abs_diff_estmt = sum(abs(t1-t2));
ABS_diff(1) = abs_diff_estmt;

t1_gauss = fit([1:17]',t1,'gauss1');
t1_bw = 2*t1_gauss.c1*sqrt(log(2));

t2_gauss = fit([1:17]',t2,'gauss1');
t2_bw = 2*t2_gauss.c1*sqrt(log(2));
bw_diff_estmt = abs(t1_bw - t2_bw);
BW_diff(1) = bw_diff_estmt;

subplot(4,2,1)
plot(1:17,t1,'k','LineWidth',3); hold on;
plot(1:17,t2,'--r','LineWidth',2); 
xlim([0 19])


% slightly shifted

t2 = [0,0,0,0.5,1,1,2,1,1,0.5,0,0,0,0,0,0,0]';

corr_estmt = corr(t1,t2);
R(2) = corr_estmt;

[~,bf1] = max(t1);
[~,bf2] = max(t2);
bf_diff_estmt = abs(bf1 - bf2);
BF_diff(2) = bf_diff_estmt;

abs_diff_estmt = sum(abs(t1-t2));
ABS_diff(2) = abs_diff_estmt;

t1_gauss = fit([1:17]',t1,'gauss1');
t1_bw = 2*t1_gauss.c1*sqrt(log(2));

t2_gauss = fit([1:17]',t2,'gauss1');
t2_bw = 2*t2_gauss.c1*sqrt(log(2));
bw_diff_estmt = abs(t1_bw - t2_bw);
BW_diff(2) = bw_diff_estmt;

subplot(4,2,3)
plot(1:17,t1,'k','LineWidth',3); hold on;
plot(1:17,t2,'--r','LineWidth',2); 
xlim([0 19])


% more shifted

t2 = [0,0,0,0,0,0.5,1,1,2,1,1,0.5,0,0,0,0,0]';

corr_estmt = corr(t1,t2);
R(3) = corr_estmt;

[~,bf1] = max(t1);
[~,bf2] = max(t2);
bf_diff_estmt = abs(bf1 - bf2);
BF_diff(3) = bf_diff_estmt;

abs_diff_estmt = sum(abs(t1-t2));
ABS_diff(3) = abs_diff_estmt;

t1_gauss = fit([1:17]',t1,'gauss1');
t1_bw = 2*t1_gauss.c1*sqrt(log(2));

t2_gauss = fit([1:17]',t2,'gauss1');
t2_bw = 2*t2_gauss.c1*sqrt(log(2));
bw_diff_estmt = abs(t1_bw - t2_bw);
BW_diff(3)

subplot(4,2,5)
plot(1:17,t1,'k','LineWidth',3); hold on;
plot(1:17,t2,'--r','LineWidth',2); 
xlim([0 19])


% most shifted

t2 = [0,0,0,0,0,0,0,0.5,1,1,2,1,1,0.5,0,0,0]';

corr_estmt = corr(t1,t2);
R(4) = corr_estmt;

[~,bf1] = max(t1);
[~,bf2] = max(t2);
bf_diff_estmt = abs(bf1 - bf2);
BF_diff(4) = bf_diff_estmt;

abs_diff_estmt = sum(abs(t1-t2));
ABS_diff(4) = abs_diff_estmt;

t1_gauss = fit([1:17]',t1,'gauss1');
t1_bw = 2*t1_gauss.c1*sqrt(log(2));

t2_gauss = fit([1:17]',t2,'gauss1');
t2_bw = 2*t2_gauss.c1*sqrt(log(2));
bw_diff_estmt = abs(t1_bw - t2_bw);
BW_diff(4) = bw_diff_estmt;

subplot(4,2,7)
plot(1:17,t1,'k','LineWidth',3); hold on;
plot(1:17,t2,'--r','LineWidth',2); 
xlim([0 19])


% plot the metrics for each of these curves

subplot(4,2,[2 4 6 8]); hold on;
plot(1:4,1-R,'LineWidth',2)
plot(1:4,BF_diff,'LineWidth',2)
plot(1:4,ABS_diff,'LineWidth',2)
plot(1:4, BW_diff,'LineWidth',2)

legend({'1-Corr','BF diff','Absolute diff','Bandwidth diff'},'Location','northwest')
set(gca,'Fontsize',16)


%% Toy example 3: differing magnitudes

figure; 

R = zeros(4,1);
BF_diff = zeros(4,1);
ABS_diff = zeros(4,1);
BW_diff = zeros(4,1);

t1 = [0,0,0,1,2,1,0,0,0,0,0,0,0,0,0,0,0]';
t2 = t1./2;

corr_estmt = corr(t1,t2);
R(1) = corr_estmt;

[~,bf1] = max(t1);
[~,bf2] = max(t2);
bf_diff_estmt = abs(bf1 - bf2);
BF_diff(1) = bf_diff_estmt;

abs_diff_estmt = sum(abs(t1-t2));
ABS_diff(1) = abs_diff_estmt;

t1_gauss = fit([1:17]',t1,'gauss1');
t1_bw = 2*t1_gauss.c1*sqrt(log(2));

t2_gauss = fit([1:17]',t2,'gauss1');
t2_bw = 2*t2_gauss.c1*sqrt(log(2));
bw_diff_estmt = abs(t1_bw - t2_bw);
BW_diff(1) = bw_diff_estmt;

subplot(4,2,1)
plot(1:17,t1,'k','LineWidth',3); hold on;
plot(1:17,t2,'--r','LineWidth',2); 
xlim([0 19])
% text(10,max(t1)-1,sprintf('R: %.2f\nBF diff: %d\nABS diff: %d\nBW diff: %0.2f',corr_estmt,bf_diff_estmt,abs_diff_estmt,bw_diff_estmt),'FontSize',18)


% slightly shifted

t2 = [0,0,0,0,0,1,2,1,0,0,0,0,0,0,0,0,0]'./2;

corr_estmt = corr(t1,t2);
R(2) = corr_estmt;

[~,bf1] = max(t1);
[~,bf2] = max(t2);
bf_diff_estmt = abs(bf1 - bf2);
BF_diff(2) = bf_diff_estmt;

abs_diff_estmt = sum(abs(t1-t2));
ABS_diff(2) = abs_diff_estmt;

t1_gauss = fit([1:17]',t1,'gauss1');
t1_bw = 2*t1_gauss.c1*sqrt(log(2));

t2_gauss = fit([1:17]',t2,'gauss1');
t2_bw = 2*t2_gauss.c1*sqrt(log(2));
bw_diff_estmt = abs(t1_bw - t2_bw);
BW_diff(2) = bw_diff_estmt;

subplot(4,2,3)
plot(1:17,t1,'k','LineWidth',3); hold on;
plot(1:17,t2,'--r','LineWidth',2); 
xlim([0 19])
% text(10,max(t1)-1,sprintf('R: %.2f\nBF diff: %d\nABS diff: %d\nBW diff: %0.2f',corr_estmt,bf_diff_estmt,abs_diff_estmt,bw_diff_estmt),'FontSize',18)



% more shifted

t2 = [0,0,0,0,0,0,0,0,1,2,1,0,0,0,0,0,0]'./2;

corr_estmt = corr(t1,t2);
R(3) = corr_estmt;

[~,bf1] = max(t1);
[~,bf2] = max(t2);
bf_diff_estmt = abs(bf1 - bf2);
BF_diff(3) = bf_diff_estmt;

abs_diff_estmt = sum(abs(t1-t2));
ABS_diff(3) = abs_diff_estmt;

t1_gauss = fit([1:17]',t1,'gauss1');
t1_bw = 2*t1_gauss.c1*sqrt(log(2));

t2_gauss = fit([1:17]',t2,'gauss1');
t2_bw = 2*t2_gauss.c1*sqrt(log(2));
bw_diff_estmt = abs(t1_bw - t2_bw);
BW_diff(3) = bw_diff_estmt;

subplot(4,2,5)
plot(1:17,t1,'k','LineWidth',3); hold on;
plot(1:17,t2,'--r','LineWidth',2); 
xlim([0 19])
% text(10,max(t1)-1,sprintf('R: %.2f\nBF diff: %d\nABS diff: %d\nBW diff: %0.2f',corr_estmt,bf_diff_estmt,abs_diff_estmt,bw_diff_estmt),'FontSize',18)


% most shifted

t2 = [0,0,0,0,0,0,0,0,0,0,0,0,1,2,1,0,0]'./2;

corr_estmt = corr(t1,t2);
R(4) = corr_estmt;

[~,bf1] = max(t1);
[~,bf2] = max(t2);
bf_diff_estmt = abs(bf1 - bf2);
BF_diff(4) = bf_diff_estmt;

abs_diff_estmt = sum(abs(t1-t2));
ABS_diff(4) = abs_diff_estmt;

t1_gauss = fit([1:17]',t1,'gauss1');
t1_bw = 2*t1_gauss.c1*sqrt(log(2));

t2_gauss = fit([1:17]',t2,'gauss1');
t2_bw = 2*t2_gauss.c1*sqrt(log(2));
bw_diff_estmt = abs(t1_bw - t2_bw);
BW_diff(4) = bw_diff_estmt;

subplot(4,2,7)
plot(1:17,t1,'k','LineWidth',3); hold on;
plot(1:17,t2,'--r','LineWidth',2); 
xlim([0 19])
% text(10,max(t1)-1,sprintf('R: %.2f\nBF diff: %d\nABS diff: %d\nBW diff: %0.2f',corr_estmt,bf_diff_estmt,abs_diff_estmt,bw_diff_estmt),'FontSize',18)



% plot the metrics for each of these curves

subplot(4,2,[2 4 6 8]); hold on;
plot(1:4,1-R,'LineWidth',2)
plot(1:4,BF_diff,'LineWidth',2)
plot(1:4,ABS_diff,'LineWidth',2)
plot(1:4, BW_diff,'LineWidth',2)

legend({'1-Corr','BF diff','Absolute diff','Bandwidth diff'},'Location','northwest')
set(gca,'Fontsize',16)


%% Toy example 4: Tuning against noise

figure; 

R = zeros(4,1);
BF_diff = zeros(4,1);
ABS_diff = zeros(4,1);
BW_diff = zeros(4,1);

t1 = [0,0,0,1,2,1,0,0,0,0,0,0,0,0,0,0,0]';
t2 = (0.05).*rand(17,1);

corr_estmt = corr(t1,t2);
R(1) = corr_estmt;

[~,bf1] = max(t1);
[~,bf2] = max(t2);
bf_diff_estmt = abs(bf1 - bf2);
BF_diff(1) = bf_diff_estmt;

abs_diff_estmt = sum(abs(t1-t2));
ABS_diff(1) = abs_diff_estmt;

t1_gauss = fit([1:17]',t1,'gauss1');
t1_bw = 2*t1_gauss.c1*sqrt(log(2));

t2_gauss = fit([1:17]',t2,'gauss1');
t2_bw = 2*t2_gauss.c1*sqrt(log(2));
bw_diff_estmt = abs(t1_bw - t2_bw);
BW_diff(1) = bw_diff_estmt;

subplot(4,2,1)
plot(1:17,t1,'k','LineWidth',3); hold on;
plot(1:17,t2,'--r','LineWidth',2); 
xlim([0 19])
% text(10,max(t1)-1,sprintf('R: %.2f\nBF diff: %d\nABS diff: %d\nBW diff: %0.2f',corr_estmt,bf_diff_estmt,abs_diff_estmt,bw_diff_estmt),'FontSize',18)


% slightly noisier

t2 = (0.2).*rand(17,1);

corr_estmt = corr(t1,t2);
R(2) = corr_estmt;

[~,bf1] = max(t1);
[~,bf2] = max(t2);
bf_diff_estmt = abs(bf1 - bf2);
BF_diff(2) = bf_diff_estmt;

abs_diff_estmt = sum(abs(t1-t2));
ABS_diff(2) = abs_diff_estmt;

t1_gauss = fit([1:17]',t1,'gauss1');
t1_bw = 2*t1_gauss.c1*sqrt(log(2));

t2_gauss = fit([1:17]',t2,'gauss1');
t2_bw = 2*t2_gauss.c1*sqrt(log(2));
bw_diff_estmt = abs(t1_bw - t2_bw);
BW_diff(2) = bw_diff_estmt;

subplot(4,2,3)
plot(1:17,t1,'k','LineWidth',3); hold on;
plot(1:17,t2,'--r','LineWidth',2); 
xlim([0 19])
% text(10,max(t1)-1,sprintf('R: %.2f\nBF diff: %d\nABS diff: %d\nBW diff: %0.2f',corr_estmt,bf_diff_estmt,abs_diff_estmt,bw_diff_estmt),'FontSize',18)



% more noisy

t2 = (0.5).*rand(17,1);

corr_estmt = corr(t1,t2);
R(3) = corr_estmt;

[~,bf1] = max(t1);
[~,bf2] = max(t2);
bf_diff_estmt = abs(bf1 - bf2);
BF_diff(3) = bf_diff_estmt;

abs_diff_estmt = sum(abs(t1-t2));
ABS_diff(3) = abs_diff_estmt;

t1_gauss = fit([1:17]',t1,'gauss1');
t1_bw = 2*t1_gauss.c1*sqrt(log(2));

t2_gauss = fit([1:17]',t2,'gauss1');
t2_bw = 2*t2_gauss.c1*sqrt(log(2));
bw_diff_estmt = abs(t1_bw - t2_bw);
BW_diff(3) = bw_diff_estmt;

subplot(4,2,5)
plot(1:17,t1,'k','LineWidth',3); hold on;
plot(1:17,t2,'--r','LineWidth',2); 
xlim([0 19])
% text(10,max(t1)-1,sprintf('R: %.2f\nBF diff: %d\nABS diff: %d\nBW diff: %0.2f',corr_estmt,bf_diff_estmt,abs_diff_estmt,bw_diff_estmt),'FontSize',18)


% noisiest

t2 = rand(17,1);

corr_estmt = corr(t1,t2);
R(4) = corr_estmt;

[~,bf1] = max(t1);
[~,bf2] = max(t2);
bf_diff_estmt = abs(bf1 - bf2);
BF_diff(4) = bf_diff_estmt;

abs_diff_estmt = sum(abs(t1-t2));
ABS_diff(4) = abs_diff_estmt;

% t1_gauss = fit([1:17]',t1,'gauss1');
% t1_bw = 2*t1_gauss.c1*sqrt(log(2));
% 
% t2_gauss = fit([1:17]',t2,'gauss1');
% t2_bw = 2*t2_gauss.c1*sqrt(log(2));
% bw_diff_estmt = abs(t1_bw - t2_bw);
% BW_diff(4) = bw_diff_estmt;

subplot(4,2,7)
plot(1:17,t1,'k','LineWidth',3); hold on;
plot(1:17,t2,'--r','LineWidth',2); 
xlim([0 19])
% text(10,max(t1)-1,sprintf('R: %.2f\nBF diff: %d\nABS diff: %d\nBW diff: %0.2f',corr_estmt,bf_diff_estmt,abs_diff_estmt,bw_diff_estmt),'FontSize',18)



% plot the metrics for each of these curves

subplot(4,2,[2 4 6 8]); hold on;
plot(1:4,1-R,'LineWidth',2)
plot(1:4,BF_diff,'LineWidth',2)
plot(1:4,ABS_diff,'LineWidth',2)
% plot(1:4, BW_diff,'LineWidth',2)

legend({'1-Corr','BF diff','Absolute diff'},'Location','northwest')
set(gca,'Fontsize',16)


%% Toy example 5: tuning curve against high-pass tuning curve

figure; 

R = zeros(4,1);
BF_diff = zeros(4,1);
ABS_diff = zeros(4,1);
BW_diff = zeros(4,1);

t1 = [0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1]'.*2;
t2 = [0,0,1,2,1,0,0,0,0,0,0,0,0,0,0,0,0]';

corr_estmt = corr(t1,t2);
R(1) = corr_estmt;

[~,bf1] = max(t1);
[~,bf2] = max(t2);
bf_diff_estmt = abs(bf1 - bf2);
BF_diff(1) = bf_diff_estmt;

abs_diff_estmt = sum(abs(t1-t2));
ABS_diff(1) = abs_diff_estmt;

t1_gauss = fit([1:17]',t1,'gauss1');
t1_bw = 2*t1_gauss.c1*sqrt(log(2));

t2_gauss = fit([1:17]',t2,'gauss1');
t2_bw = 2*t2_gauss.c1*sqrt(log(2));
bw_diff_estmt = abs(t1_bw - t2_bw);
BW_diff(1) = bw_diff_estmt;

subplot(4,2,7)
plot(1:17,t1,'k','LineWidth',3); hold on;
plot(1:17,t2,'--r','LineWidth',2); 
xlim([0 19])
% text(10,max(t1)-1,sprintf('R: %.2f\nBF diff: %d\nABS diff: %d\nBW diff: %0.2f',corr_estmt,bf_diff_estmt,abs_diff_estmt,bw_diff_estmt),'FontSize',18)


% slightly shifted

t2 = [0,0,0,0,0,1,2,1,0,0,0,0,0,0,0,0,0]';

corr_estmt = corr(t1,t2);
R(2) = corr_estmt;

[~,bf1] = max(t1);
[~,bf2] = max(t2);
bf_diff_estmt = abs(bf1 - bf2);
BF_diff(2) = bf_diff_estmt;

abs_diff_estmt = sum(abs(t1-t2));
ABS_diff(2) = abs_diff_estmt;

t1_gauss = fit([1:17]',t1,'gauss1');
t1_bw = 2*t1_gauss.c1*sqrt(log(2));

t2_gauss = fit([1:17]',t2,'gauss1');
t2_bw = 2*t2_gauss.c1*sqrt(log(2));
bw_diff_estmt = abs(t1_bw - t2_bw);
BW_diff(2) = bw_diff_estmt;

subplot(4,2,5)
plot(1:17,t1,'k','LineWidth',3); hold on;
plot(1:17,t2,'--r','LineWidth',2); 
xlim([0 19])
% text(10,max(t1)-1,sprintf('R: %.2f\nBF diff: %d\nABS diff: %d\nBW diff: %0.2f',corr_estmt,bf_diff_estmt,abs_diff_estmt,bw_diff_estmt),'FontSize',18)



% more shifted

t2 = [0,0,0,0,0,0,0,0,1,2,1,0,0,0,0,0,0]';

corr_estmt = corr(t1,t2);
R(3) = corr_estmt;

[~,bf1] = max(t1);
[~,bf2] = max(t2);
bf_diff_estmt = abs(bf1 - bf2);
BF_diff(3) = bf_diff_estmt;

abs_diff_estmt = sum(abs(t1-t2));
ABS_diff(3) = abs_diff_estmt;

t1_gauss = fit([1:17]',t1,'gauss1');
t1_bw = 2*t1_gauss.c1*sqrt(log(2));

t2_gauss = fit([1:17]',t2,'gauss1');
t2_bw = 2*t2_gauss.c1*sqrt(log(2));
bw_diff_estmt = abs(t1_bw - t2_bw);
BW_diff(3) = bw_diff_estmt;

subplot(4,2,3)
plot(1:17,t1,'k','LineWidth',3); hold on;
plot(1:17,t2,'--r','LineWidth',2); 
xlim([0 19])
% text(10,max(t1)-1,sprintf('R: %.2f\nBF diff: %d\nABS diff: %d\nBW diff: %0.2f',corr_estmt,bf_diff_estmt,abs_diff_estmt,bw_diff_estmt),'FontSize',18)


% most shifted

t2 = [0,0,0,0,0,0,0,0,0,0,0,0,1,2,1,0,0]';

corr_estmt = corr(t1,t2);
R(4) = corr_estmt;

[~,bf1] = max(t1);
[~,bf2] = max(t2);
bf_diff_estmt = abs(bf1 - bf2);
BF_diff(4) = bf_diff_estmt;

abs_diff_estmt = sum(abs(t1-t2));
ABS_diff(4) = abs_diff_estmt;

t1_gauss = fit([1:17]',t1,'gauss1');
t1_bw = 2*t1_gauss.c1*sqrt(log(2));

t2_gauss = fit([1:17]',t2,'gauss1');
t2_bw = 2*t2_gauss.c1*sqrt(log(2));
bw_diff_estmt = abs(t1_bw - t2_bw);
BW_diff(4) = bw_diff_estmt;

subplot(4,2,1)
plot(1:17,t1,'k','LineWidth',3); hold on;
plot(1:17,t2,'--r','LineWidth',2); 
xlim([0 19])
% text(10,max(t1)-1,sprintf('R: %.2f\nBF diff: %d\nABS diff: %d\nBW diff: %0.2f',corr_estmt,bf_diff_estmt,abs_diff_estmt,bw_diff_estmt),'FontSize',18)



% plot the metrics for each of these curves

subplot(4,2,[2 4 6 8]); hold on;
plot(1:4,1-R,'LineWidth',2)
plot(1:4,BF_diff,'LineWidth',2)
plot(1:4,ABS_diff,'LineWidth',2)
plot(1:4, BW_diff,'LineWidth',2)

legend({'1-Corr','BF diff','Absolute diff','Bandwidth diff'},'Location','northwest')
set(gca,'Fontsize',16)
