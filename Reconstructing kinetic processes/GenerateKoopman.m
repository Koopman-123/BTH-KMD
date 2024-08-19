
function [Psi] = GenerateKoopman(data, specified_modes)

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
Avg = mean(Data, 2); % Compute time average
[eigval, Modes1, bo] = H_DMD(Data - repmat(Avg, 1, size(Data, 2)), delay);
toc;

disp('Sorting Modes...');
tic;
omega = log(diag(eigval)) / delt; % Continuous time eigenvalues
decayRates = -real(omega); % Decay rates
growthRates = real(omega); % Growth rates
[T, Im] = sort(abs(omega), 'descend'); % Sort modes by the magnitude of growth/decay rates
Modes1 = Modes1(:, Im);
bo = bo(Im);
toc;
    %% neutral+unstable+stable
%     neutral_eig = eigval(abs((abs(eigval)-1))<0.001);neutral_idx = find(abs((abs(eigval)-1))<0.001);
%     unstable_eig = eigval((abs(eigval)-1)>0.001);stable_idx = find((abs(eigval)-1)>0.001);
%     stable_eig = eigval((abs(eigval)-1)<-0.001);unstable_idx = find((abs(eigval)-1)<-0.001);
%     
%     aa = real(diag(neutral_eig)); % Real part of eigenvalues
%     bb = imag(diag(neutral_eig)); % Imaginary part of eigenvalues
%     
%     figure()
%     scatter(real(neutral_eig), imag(neutral_eig), 10, 'b', 'filled');
%     hold on
%     scatter(real(stable_eig), imag(stable_eig), 10, 'k',  'filled');
%     hold on 
%     scatter(real(unstable_eig), imag(unstable_eig), 10,'r',  'filled');
%     legend('Netural','Stable','Unstable')
%     hold off 
%% Compute and Plot Specified Modes
disp('Computing and Plotting Specified Modes...')
tic;
[nbx, nbt] = size(Data); % Get data size
time = (0:nbt - 1) * delt; % Specify time interval
Psi = zeros(nbx, nbt, length(specified_modes));

X_reconstructed = zeros(nbx, nbt); % Initialize the reconstructed data matrix
mode_energies = zeros(1, length(specified_modes)); 
folder = 'D:\\Desktop\\BTH-PM2.5\\Reconstructing kinetic processes\\';

for idx = 1:length(specified_modes)
    mode_index = specified_modes(idx);
    psi = zeros(1, nbt); % Preallocate time evolution of mode
    omeganow = omega(mode_index); % Get current eigenvalue
    bnow = bo(mode_index); % Get current amplitude coefficient
    for t = 1:length(time)
        psi(:, t) = exp(omeganow * time(t)) * bnow; % Compute mode's time evolution
    end
    psi = Modes1(1:nbx, mode_index) * psi;
    Psi(:, :, idx) = psi;
    
    X_reconstructed = X_reconstructed + psi; % Accumulate the mode contributions
    mode_energies(idx) = norm(psi, 2)^2;
end
    
mmm = abs(psi);   
Avg_matrix = repmat(Avg, 1, nbt); % Extend Avg to match the dimensions of X_reconstructed
X_reconstructed = real(X_reconstructed) + Avg_matrix; 
csvwrite([folder, 'X_reconstructed_stable.csv'], X_reconstructed);
csvwrite([folder, 'X_reconstructed_unstable.csv'], X_reconstructed);
csvwrite([folder, 'X_reconstructed_neutral.csv'], X_reconstructed);
csvwrite([folder, 'X_reconstructed_stable_unstable.csv'], X_reconstructed);
csvwrite([folder, 'X_reconstructed_stable_unstable_neutral.csv'], X_reconstructed);
%     csvwrite([folder, 'mode_energies.csv'], mode_energies'); % Save mode energies
%     csvwrite([folder, 'growth_rates.csv'], decayRates);
toc
end