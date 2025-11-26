% Plot_Density_2D_Fourier.m
% Clear workspace and command window
clear;
clc;

% -------------------------------
% Parameters (from Config)
% -------------------------------
cfg = GetCommonConfig();

% Use the first D value if multiple are provided
Db = cfg.D_values(1);

% -------------------------------
% Determine filename based on filling type
% -------------------------------
FileNamePostfix = GetSimulationParamsPostfix(cfg.Ly, cfg.Lx, cfg.t2, cfg.J, ...
    cfg.FillingType, cfg.FillingValue, Db, cfg.Geometry);

% Full file path (adjust if needed)
charge_density_file = ['../data/nf', FileNamePostfix];

% -------------------------------
% Load charge density data from JSON
% -------------------------------
try
    ChargeDensityData = jsondecode(fileread(char(charge_density_file)));
catch ME
    error('Failed to load charge density data (%s): %s', charge_density_file, ME.message);
end

% -------------------------------
% Extract coordinates and densities
% -------------------------------
N = size(ChargeDensityData, 1);  % Number of data points

x_coor = zeros(1, N);
y_coor = zeros(1, N);
ChargeDensity = zeros(1, N);

if iscell(ChargeDensityData)
    for i = 1:N
        idx = ChargeDensityData{i}{1};
        x_coor(i) = fix(idx / cfg.Ly);
        y_coor(i) = mod(idx, cfg.Ly);
        ChargeDensity(i) = ChargeDensityData{i}{2}(1);
    end
else
    idx = ChargeDensityData(:, 1);
    x_coor = fix(idx / cfg.Ly);
    y_coor = mod(idx, cfg.Ly);
    ChargeDensity = ChargeDensityData(:, 2);
end

% -------------------------------
% Compute hole density and reshape into 2D lattice
% -------------------------------
hole_density = 1 - ChargeDensity;
hole_matrix = reshape(hole_density, [cfg.Ly, cfg.Lx]);  % both DMRG and matlab uses column-major

% Optional: Subtract mean to remove DC component (k=0 peak)
% hole_matrix = hole_matrix - mean(mean(hole_matrix));

hole_alongx = mean(hole_matrix);
% hole1D_fft = fftshift(fft(hole_alongx)); % Unused in plot currently

% -------------------------------
% Compute 2D Fourier Transform and shift zero frequency to center
% -------------------------------
hole_fft = fftshift(fft2(hole_matrix));  % 2D FFT and shift zero frequency to center
S_k = abs(hole_fft).^2;  % Power spectrum (structure factor)

% -------------------------------
% Define k-space grid for visualization
% -------------------------------
dkx = 2 * pi / cfg.Lx;
dky = 2 * pi / cfg.Ly;

kx = (-pi : dkx : pi - dkx);  % From -π to π 
ky = (-pi : dky : pi - dky);

% Ensure grid matches S_k size exactly (handling odd/even sizes roughly)
% A more robust way if sizes don't match is to linspace based on size(S_k)
if length(kx) ~= size(S_k, 2) || length(ky) ~= size(S_k, 1)
   kx = linspace(-pi, pi - 2*pi/cfg.Lx, size(S_k, 2));
   ky = linspace(-pi, pi - 2*pi/cfg.Ly, size(S_k, 1));
end

[KX, KY] = meshgrid(kx, ky);  % Create meshgrid for k-space

% -------------------------------
% Plotting the structure factor in the First Brillouin Zone
% -------------------------------
figure;
pcolor(KX, KY, S_k);
shading interp;
colorbar;
title('Structure Factor of Hole Density in First Brillouin Zone');
xlabel('k_x (\pi)');
ylabel('k_y (\pi)');
xlim([-pi pi]);
ylim([-pi pi]);
axis equal;
colormap('hot');  % Use 'hot' for better contrast

% Optional: Log scale for enhanced contrast
% S_k_log = log(1 + S_k);  % Avoid log(0)
% figure;
% pcolor(KX, KY, S_k_log);
% shading interp;
% colorbar;
% title('Structure Factor of Hole Density (Log Scale)');
% xlabel('k_x (\pi)');
% ylabel('k_y (\pi)');
% xlim([-pi pi]);
% ylim([-pi pi]);
% axis equal;
% colormap('hot');

hold off;
