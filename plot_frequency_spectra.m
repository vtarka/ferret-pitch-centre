% AUTHOR: Veronica Tarka, veronica.tarka@dpag.ox.ac.uk, May 2023

figure;
hold on
set(gca,'fontsize',22)
xlabel('Pitch (Hz)')
ylabel('Frequency Spectrum (Hz)')

CT = {500:500:28000, 1000:1000:28000, 2000:2000:28000, 4000:4000:28000 };
low = {1000:500:2000, 2000:1000:6000, 4000:2000:14000,8000:4000:16000};
high = {3000:500:28000, 8000:1000:28000, 18000:2000:28000,[24000 28000]};

pitches = [500 1000 2000 4000];

width = 0.7;
gap = 0.1;
big_gap = 1.3;

line_widths = cell(length(pitches),1);
xtks = zeros(length(pitches),1);

starting_point = 0.5;
for pp = 1:length(pitches)
    lines = zeros(3,2);

    for i = 1:3
        lines(i,1) = (width*(i-1)) + gap*(i-1) + starting_point;
        lines(i,2) = (width*i) + gap*(i-1) + starting_point;
    end

    line_widths{pp} = lines;

    xtks(pp) = (lines(1,1) + lines(3,2)) / 2;
    starting_point = lines(3,2) + big_gap;
end

xlim([0 starting_point])
xticks(xtks)
xticklabels(pitches)


for pp = 1:length(pitches)
    line_width = line_widths{pp};
    line_width = line_width(1,:);

    freqs = CT{pp};

    for ff = 1:length(freqs)
        line(line_width,[freqs(ff),freqs(ff)],'Color','black','linewidth',2)
    end
end

xticks()
pause

for pp = 1:length(pitches)

    line_width = line_widths{pp};
    line_width = line_width(2,:);

    freqs = low{pp};

    for ff = 1:length(freqs)
        line(line_width,[freqs(ff),freqs(ff)],'Color','blue','linewidth',3)
    end
end

pause

for pp = 1:length(pitches)
    
    line_width = line_widths{pp};
    line_width = line_width(3,:);

    freqs = high{pp};

    for ff = 1:length(freqs)
        line(line_width,[freqs(ff),freqs(ff)],'Color','red','linewidth',3)
    end
end