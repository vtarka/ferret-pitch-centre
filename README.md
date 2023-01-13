# ferret-pitch-centre

# Useful info on data formatting
Each of spiking file contains 5 variables:

    Y - n x 6 matrix, column 2 is the spike time relative to stimulus onset, column 3 is which unit was spiking, column 4 is the stimulus number (name for the unique combo of sound type (CT0, pure tone, etc) and the F0), column 5 denotes the repeat of the stimulus number, and column 6 is the trial number.
    type - 217 x 1 vector, labels the stimulus type corresponding to the stimulus number in column 4 of Y. For example, if column 4 holds the number 198, that spike corresponds to a trial where type(198) was played.
    F0 - 217 x 1 vector, labels the F0 of the stimulus number in column 4 of Y. Same idea as the type variable... If column 4 of Y holds 198, that spike corresponds to a trial where F0(198) Hz was the fundamental.
    sensitivity - nUnits x 13 matrix, each index labels whether the unit was determined to be frequency sensitive to that particular stimulus type (1 is yes), sensitivity(1,:) is the labelling for the first unit in unique(Y(:,3)) and so on...
    BFs - nUnits x 13 matrix, each index labels which F0 evoked the maximum firing rate in the averaged tuning curve for that particular stimulus type, an entry of 0 means the unit wasn't frequency-sensitive to that stimulus type. If BFs(1,13) = 7, that means the 1st unit in unique(Y(:,3)) had maximum firing at the 7th entry in unique(F0) in response to pure tones.

The stimulus types are labeled as:
1 - CT0
2 - CT10
3 - CT20
4 - CT40
5 - CT5
6 - F0MaskHigh
7 - F0MaskLow
8 - allHarm
9 - alt
10 - high
11 - low
12 - rand
13 - tone
