% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, May 2023


parent_dir = '/media/veronica/Kat Data/Veronica/pitch_ephys/Quentin Pitch ephys/stim2020/70/';
folder = '/media/veronica/Kat Data/Veronica/pitch_ephys/Quentin Pitch ephys/stim2020/70/*.wav';

listing = dir(folder);

for ff = 1:length(listing)

    [data,fs] = audioread([parent_dir listing(ff).name]);
    T = 1/fs;
    L = length(data);

    data_fft = fft(data);
    P2 = abs(data_fft/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(L/2))/L;

    plot(f,P1)
    title(listing(ff).name)

    pause
end