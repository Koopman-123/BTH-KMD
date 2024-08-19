function [Psi] = GenerateKoopmanModes(data, mode1, mode2, month)
    %% Load Data
    disp('Loading Data Set...');
    tic;
    if strcmp(data, 'Day_mean');
    Data = dlmread('daily PM2.5.txt');
    delay = 22; dtype = 'Mean'; delt = 1; delx = 1;
    hwy = 'day'; hwylength = 76; xpath = 'x76.txt'; ypath = 'y76.txt';
    end
    toc;
    %% Compute KMD and Sort Modes 
    disp('Computing KMD via Hankel-DMD...');
    tic;
    Avg = mean(Data, 2);
    [eigval, Modes1, bo] = H_DMD(Data - repmat(Avg, 1, size(Data, 2)), delay);
    toc;
    disp('Sorting Modes...')
    tic;
    
    omega = log(diag(eigval)) / delt; % Continuous time eigenvalues
    decayRates = -real(omega); % Decay rates
    growthRates = real(omega); % Growth rates
    Freal = imag(omega) / (2 * pi); % Frequencies
    [T, Im] = sort((1 ./ Freal), 'descend'); % Sort frequencies
    omega = omega(Im); Modes1 = Modes1(:, Im); bo = bo(Im); % Sort modes
    toc;
    %% Compute and Plot Modes
    disp('Computing and Plotting Modes...')
    tic;
    [nbx, nbt] = size(Data); % Get data size
    time = (0:nbt - 1) * delt; % Specify time interval
    Psi = zeros(nbx, nbt, mode2 - mode1 + 1);
    
    % Compute amplitude coefficients for each mode
%     amplitudes = abs(bo); % Assuming 'bo' contains the amplitude coefficients
%     res = [];
    X_reconstructed = zeros(nbx, nbt); % Initialize the reconstructed data matrix
    folder = 'D:\\Desktop\\BTH-PM2.5\\Reconstructing kinetic processes\\';
%     csvwrite([folder, 'neutral.csv'], neutral_eig);
%     csvwrite([folder, 'unstable.csv'], unstable_eig);
%     csvwrite([folder, 'stable.csv'], stable_eig);
    
    for i = mode1:mode2; % Loop through all modes to plot
        psi = zeros(1, nbt); % Preallocate time evolution of mode
        omeganow = omega(i); % Get current eigenvalue 
        bnow = bo(i); % Get current amplitude coefficient 
        parfor t = 1:length(time);  
            psi(:, t) = exp(omeganow * time(t)) * bnow; 
        end
        psi = Modes1(1:nbx, i) * psi; % Compute mode 
        Psi(:, :, i) = psi; % Store & output modes
        X_reconstructed = X_reconstructed + psi; % Accumulate the mode contributions
%         amplitudes = abs(psi(:));
%         sortmode = sort(amplitudes, 'descend'); % Calculate and store mode energy
%         mode_energies(i - mode1 + 1) = sum(amplitudes.^2); % Calculate and store mode energy
%         X_reconstructed(X_reconstructed < 0) = 0; % Ensure all values are non-negative after each addition
%         X_reconstructed_abs = abs(X_reconstructed);  
%         csvwrite([folder, 'mode_', num2str(i), '.csv'], mmm); % Save mode to CSV
       end
%% Load Mode Energies from CSV
folder = 'D:\Desktop\\BTH-PM2.5\Reconstructing kinetic processes\';
% filename = [folder, 'mode_energies.csv'];
% mode_energies = csvread(filename); 
%% Calculate Total Energy
% total_energy = sum(mode_energies);
%% Calculate Cumulative Energies
% cumulative_energies = cumsum(mode_energies);
%% Calculate Cumulative Energy Percentages
% cumulative_energy = (cumulative_energies / total_energy) * 100;
%% Optionally, save the cumulative energy percentages to a CSV file
aaaa = real(psi)
bbbb = abs(psi)
Avg_matrix = repmat(Avg, 1, nbt); % Extend Avg to match the dimensions of X_reconstructed
X_reconstructed = real(X_reconstructed) + Avg_matrix; 
csvwrite([folder, 'X_reconstructed+key+mode.csv'], X_reconstructed);
% csvwrite([folder, 'growth_rates.csv'], decayRates);
% csvwrite([folder, 'mode_energies.csv'], mode_energies); 
% csvwrite([folder, 'energy_percentages.csv'], energy_percentages); 
% csvwrite([folder, 'cumulative_energy_percentages.csv'], cumulative_energy_percentages); % Save cumulative energy percentages
toc
end