% Plot_Density_2D_color_LocalSpin.m
% ====== Model Parameters from Config ===== %
clear;
cfg = GetCommonConfig();

% Visualization Parameters (Local)
max_arrow_length = 0.4;

% Generate Filename
% Use first D value if multiple are present
Db = cfg.D_values(1);

FileNamePostfix = GetSimulationParamsPostfix(cfg.Ly, cfg.Lx, cfg.t2, cfg.J, ...
    cfg.FillingType, cfg.FillingValue, Db, cfg.Geometry);

try
    ChargeDensityData = jsondecode(fileread(['../data/nf', FileNamePostfix]));
    SpinData = jsondecode(fileread(['../data/sz', FileNamePostfix]));
catch ME
    error('Failed to load data (Postfix: %s): %s', FileNamePostfix, ME.message);
end

% Extract coordinates and densities (unchanged logic)
x_coor = zeros(1, size(ChargeDensityData, 1));
y_coor = zeros(1, size(ChargeDensityData, 1));
ChargeDensity = zeros(1, numel(x_coor));
SpinDensity = zeros(1, numel(x_coor));

if iscell(ChargeDensityData)
    for i = 1:numel(x_coor)
        x_coor(i) = fix(ChargeDensityData{i}{1} / cfg.Ly);
        y_coor(i) = mod(ChargeDensityData{i}{1}, cfg.Ly);
        ChargeDensity(i) = transpose(ChargeDensityData{i}{2}(1));
        SpinDensity(i) = transpose(SpinData{i}{2}(1));
    end
else
    x_coor = fix(ChargeDensityData(:,1) / cfg.Ly);
    y_coor = mod(ChargeDensityData(:,1), cfg.Ly);
    ChargeDensity = (ChargeDensityData(:,2));
    SpinDensity = (SpinData(:,2));
end

% Calculate hole density and reshape matrices
hole_density = 1 - ChargeDensity;
ave_hole_density = mean(hole_density);
fprintf("ave hole density : %.12f\n", ave_hole_density);
hole_matrix = reshape(hole_density, [cfg.Ly, cfg.Lx]);
spin_matrix = reshape(SpinDensity, [cfg.Ly, cfg.Lx]);

% Create figure and plot colormap
figure;
imagesc(1:cfg.Lx, 1:cfg.Ly, hole_matrix);
set(gca, 'YDir', 'reverse');
try
    colormap(slanCM(19)); 
catch
    colormap('parula'); % Fallback if slanCM is missing
end
colorbar;
hold on;

% Calculate scaled spin arrows
max_spin_abs = max(abs(spin_matrix(:)));
if max_spin_abs > 0
    V_scaled = spin_matrix * (max_arrow_length / max_spin_abs);
else
    V_scaled = zeros(size(spin_matrix));
end

% Generate grid for arrows
[X, Y] = meshgrid(1:cfg.Lx, 1:cfg.Ly);

% Plot positive spins (downward, blue)
positive = spin_matrix > 0;
v_pos = -V_scaled(positive);
quiver(X(positive), Y(positive) - v_pos/2, zeros(nnz(positive), 1), v_pos, ...
    0, 'Color', 'b', 'MaxHeadSize', 1, 'LineWidth', 1.5, 'AutoScale', 'off');

% Plot negative spins (upward, red)
negative = spin_matrix < 0;
v_neg = -V_scaled(negative);
quiver(X(negative), Y(negative) - v_neg/2, zeros(nnz(negative), 1), v_neg, ...
    0, 'Color', 'r', 'MaxHeadSize', 1, 'LineWidth', 1.5, 'AutoScale', 'off');


spin_flat = spin_matrix(:);
digits = 1;   % Number of significant figures to keep
max_spin = max(abs(spin_flat));
if max_spin > 0
    demonstrate_spin = round(max_spin, digits - floor(log10(abs(max_spin))) - 1);
else
    demonstrate_spin = 0;
end

% Spin legend
legend_x = cfg.Lx / 2;  % Position to the right of the plot
legend_y_start = cfg.Ly + 2;  % Start position in the middle of the plot
legend_spacing = 1;  % Vertical spacing between legend items

% Calculate the scaled arrow length for the legend
if max(abs(spin_flat)) > 0
    legend_arrow_length = max_arrow_length * demonstrate_spin / max(abs(spin_flat));
else
    legend_arrow_length = 0;
end

% Spin Up example
v_legend_up = -legend_arrow_length;
quiver(legend_x, legend_y_start - v_legend_up/2, ...
    0, v_legend_up, ...
    0, 'Color', 'b', 'MaxHeadSize', 1, 'LineWidth', 1.5);
text(legend_x + 0.5, legend_y_start, ...
    ['Spin Up = ', num2str(demonstrate_spin)], 'FontSize', 14, 'VerticalAlignment', 'middle');

% Spin Down example
v_legend_down = legend_arrow_length;
quiver(legend_x, legend_y_start - legend_spacing - v_legend_down/2, ...
    0, v_legend_down, ...
    0, 'Color', 'r', 'MaxHeadSize', 1, 'LineWidth', 1.5);
text(legend_x + 0.5, legend_y_start - legend_spacing, ...
    ['Spin Down = ', num2str(demonstrate_spin)], 'FontSize', 14, 'VerticalAlignment', 'middle');

hold off;

% Adjust axis properties
axis equal;
xlim([0.5, cfg.Lx + 0.5]);
ylim([0.5, cfg.Ly + 3]);
set(gca, 'YDir', 'reverse');
axis off;

% Label colorbar
c = colorbar;
ylabel(c, 'Hole Density', 'FontSize',14);
