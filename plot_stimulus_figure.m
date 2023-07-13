% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, May 2023

stims = {'low','high','alt','rand','allHarm','CT0','CT5','CT10','CT20','CT40'};

sps = [1 6; 2 7; 3 8; 4 9; 5 10; 11 16; 12 17; 13 18; 14 19; 15 20];

figure;
color = 'k';

for ss = 1:length(stims)

%     if strcmp(color,'k')
%         color = 'k';
%     else
%         color = 'k';
%     end

    if ss == 1
        color = 'b';
    elseif ss == 2
        color = 'r';
    else
        color = 'k';
    end

    [data,fs] = audioread(['/media/veronica/Kat Data/Veronica/pitch_ephys/Quentin Pitch ephys/stim2020/70/Pitch2018_' stims{ss} '_2000Hz_70dB.wav']);

    T = 1/fs;
    L = length(data);

    data_fft = fft(data);
    P2 = abs(data_fft/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(L/2))/L;

    P1(f>30000) = [];
    f(f>30000) = [];

    subplot(4,5,sps(ss,1));
    plot(f,P1,color,'linewidth',2)
    xticks(0:10000:30000)
    xticklabels(0:10:30)
    yticks([max(P1)/2 max(P1)])
    yticklabels([-75 -50])
    ylim([0 max(P1)+0.0001])
    set(gca,'fontsize',18)

    stream_sample = data(10000:10500);
    subplot(4,5,sps(ss,2));
    plot(stream_sample,color)
    xlim([0 500])
    xticks(0:250:500)
    xticklabels(0:2.5:5)
    yticks([min(stream_sample) 0 max(stream_sample)])
    yticklabels([-1 0 1])
    ylim([min(stream_sample)-0.005 max(stream_sample)+0.005])
    set(gca,'fontsize',18)

end