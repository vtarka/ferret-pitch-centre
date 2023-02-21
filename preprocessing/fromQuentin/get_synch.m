function synch_ch = get_synch(data_root) 
%% Load the synch channel
fprintf('== Loading the synch channel ==\n');tic;

dir_info = dir([data_root '/synch*.mat']);

if ~isempty(dir_info)
    disp('syn channel already extracted.')
    load(fullfile(dir_info.folder,dir_info.name));
    synch_ch = synch;
    
else
    % Get syn channel from data file
    dir_info = dir([data_root '/*ap.bin']); %Get the names of all files with this extension
    % bin_filename = fullfile(dir_info.folder, dir_info.name); %Form the name of the file to be opened

    bin_filename = fullfile(data_root,dir_info.name);
    nChansTotal = 385;
    bytes_sample = 2;
    nSampsTotal = dir_info.bytes/nChansTotal/bytes_sample; %Find the total number of samples in the file to be analyzed

    offset_synch = 384*bytes_sample; %This should move the pointer to the 385 channel i.e. the synch channel

    fid = fopen(bin_filename, 'r'); %Open an fid for the file of interest. This is something like a pointer for this file

    fseek(fid, offset_synch, 'bof'); %Move the place to which fseek is pointing

    precision = '*int16';
    skip = 384*bytes_sample;
    synch_ch = fread(fid,[nSampsTotal,1],precision,skip);
    % synch_ch = fread(fid,[600000,1],precision,skip);
    fclose(fid); %Close the fid so it doesn't interfere with future calls of fopen
end
fprintf('== Done! Loading took %.0fs ==\n',toc);