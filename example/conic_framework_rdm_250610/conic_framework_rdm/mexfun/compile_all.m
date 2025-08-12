% Get a list of all .c files in the directory
cFiles = dir(fullfile('.', '*.c'));

% Compile each .c file with -largeArrayDims
for k = 1:length(cFiles)
    % Full path to the .c file
    cFilePath = fullfile(cFiles(k).folder, cFiles(k).name);
    
    % Compile the .c file with -largeArrayDims
    fprintf('Compiling %s...\n', cFiles(k).name);
    mex('-R2018a', cFilePath);
end

disp('All files compiled successfully.');
