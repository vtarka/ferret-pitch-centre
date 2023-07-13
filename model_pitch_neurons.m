% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, May 2023

figure; hold on

peaks = 200:250:20000;

for pp = 1:length(peaks)

    model_peak = peaks(pp);
    model_width = 10;
    
    max_freq = 30000;
    min_freq = 0;
    
    frequencies = [354 420 500 595 707 841 1000 1189 1414 1682 2000 2378 2828 3364 4000];
    
    CT_tuning_curve = [];
    low_tuning_curve = [];
    high_tuning_curve = [];
    
    for ff = 1:length(frequencies)
    
        freq = num2str(frequencies(ff));
    
        [CT_stream,fs] = audioread(['/media/veronica/Kat Data/Veronica/pitch_ephys/Quentin Pitch ephys/stim2020/70/Pitch2018_CT0_' freq 'Hz_70dB.wav']);
        low_stream = audioread(['/media/veronica/Kat Data/Veronica/pitch_ephys/Quentin Pitch ephys/stim2020/70/Pitch2018_low_' freq 'Hz_70dB.wav']);
        high_stream = audioread(['/media/veronica/Kat Data/Veronica/pitch_ephys/Quentin Pitch ephys/stim2020/70/Pitch2018_high_' freq 'Hz_70dB.wav']);
    
        T = 1/fs;
        L = length(CT_stream);
    
        CT_fft = fft(CT_stream);
        P2 = abs(CT_fft/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
    %     P1(round(max_freq/fs):end) = [];
        CT_fft = P1;
    
        low_fft = fft(low_stream);
        P2 = abs(low_fft/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
    %     P1(round(max_freq/fs):end) = [];
        low_fft = P1;
    
        high_fft = fft(high_stream);
        P2 = abs(high_fft/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
    %     P1(round(max_freq/fs):end) = [];
        high_fft = P1;
       
        f = fs*(0:(L/2))/L;
    
        tuning_model = exp(-((f-model_peak)/model_width).^2/2);
    
        CT_tuning = sum(CT_fft.*tuning_model');
        low_tuning = sum(low_fft.*tuning_model');
        high_tuning = sum(high_fft.*tuning_model');
    
        CT_tuning_curve(end+1) = CT_tuning;
        low_tuning_curve(end+1) = low_tuning;
        high_tuning_curve(end+1) = high_tuning;
    
    %     plot(ff,CT_tuning,'k.','markersize',20);
    %     plot(ff,low_tuning,'r.','markersize',20);
    %     plot(ff,high_tuning,'b.','markersize',20);
    
    end
    
    clf; hold on

    plot(CT_tuning_curve,'k','linewidth',2)
    plot(low_tuning_curve,'b','linewidth',2)
    plot(high_tuning_curve,'r','linewidth',2)
    xticks(2:2:length(frequencies))
    xticklabels(frequencies(2:2:length(frequencies)))

    title(model_peak)
    pause

end
