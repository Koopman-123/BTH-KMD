% Dynamic Mode Decomposition of PM2.5

%% Instructions for Koopman Modes.

%Plots Generated:
% The funciton will generate plots of the desired Koopman Modes.

% Examples:neutral+unstable+stable
filePath = 'D:\Desktop\BTH-PM2.5\Reconstructing kinetic processes\stable_data.csv';
dataTable = readtable(filePath);
specified_modes = dataTable.mode;
neutral_unstable_stable = GenerateKoopman('Day_mean', specified_modes);

key_mode=GenerateKoopmanModes('Day_mean', 3, 3);





