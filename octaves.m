
%% Calculates the octave difference between two frequencies
% AUTHOR: Kerry Walker

function n=octaves(f1,f2)
%calculates the number of octaves a frequency range corresponds to.

n=(log10(f1/f2))/(log10(2));
n=abs(n);
end