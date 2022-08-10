% Extracts and writes into separate files data for harmonics Z0, K1, O1, S2, M2.
% Processes all the files in the specified directory.
% There must be no other files in this directory.
% m-file: extract_exam.m
% M.Krassovski& peter, January/16/2004
% krassovskim@dfo-mpo.gc.ca

headerLines = 4;
files = dir('E:\matlab7\work\jeju6403');
numFiles = length(files);
% Remove two psedofiles (. and ..) from the file list
files = files(3:numFiles);
numFiles = numFiles - 2;
for i = 1:numFiles
%    files(i).name   %%%%%%%%%%% show current file name
    fid = fopen(['E:\matlab7\work\jeju6403', '\', files(i).name]);
    fout = fopen(['E:\matlab7\work\jeju6403', '\', 'out_', files(i).name], 'w');
    fprintf(fout, '%s\r\n', '         A      G');    % write header
    
    for skip = 1:headerLines      % skip file header
        line = fgetl(fid);
    end
    while 1
        line = fgetl(fid);              % read next line
        if line == -1, break, end       % until the end of file

        if strcmp(line(6:7), 'Z0') | strcmp(line(6:7), 'K1') | strcmp(line(6:7), 'O1') | strcmp(line(6:7), 'S2') | strcmp(line(6:7), 'M2')
            fprintf(fout, '%s\r\n', [line(6:9), line(41:53)]);
        end
    end
    fclose(fid);
    fclose(fout);
end