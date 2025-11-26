% Plot_Density_2D_Scatter_LocalSpin.m
clear;
close all;

% ====== Configuration ======
cfg = GetCommonConfig();

% Visualization Parameters
max_size = 200;      % Maximum circle size for hole density
max_arrow_length = 0.4; 
digits = 1;          % Significant figures for legend

% Use the first D value if multiple are provided
Db = cfg.D_values(1);

% Generate Filename Postfix
FileNamePostfix = GetSimulationParamsPostfix(cfg.Ly, cfg.Lx, cfg.t2, cfg.J, ...
    cfg.FillingType, cfg.FillingValue, Db, cfg.Geometry);

% ========== Load Data ==========
try
    % Load charge density
    charge_density_file = ['../data/nf', FileNamePostfix];
    ChargeDensityData = jsondecode(fileread(charge_density_file));
    
    % Load spin data
    spin_file = ['../data/sz', FileNamePostfix];
    SpinData = jsondecode(fileread(spin_file));
catch ME
    error('Failed to load data (Postfix: %s): %s', FileNamePostfix, ME.message);
end

% Extract coordinates and densities
num_sites = size(ChargeDensityData, 1);
x_coor = zeros(1, num_sites);
y_coor = zeros(1, num_sites);
ChargeDensity = zeros(1, num_sites);
SpinDensity = zeros(1, num_sites);

if iscell(ChargeDensityData) % complex/cell format
    for i = 1:num_sites
        % Coordinates
        idx = ChargeDensityData{i}{1};
        x_coor(i) = fix(idx / cfg.Ly);
        y_coor(i) = mod(idx, cfg.Ly);
        
        % Values
        ChargeDensity(i) = transpose(ChargeDensityData{i}{2}(1));
        SpinDensity(i) = transpose(SpinData{i}{2}(1));
    end
else % matrix format
    x_coor = fix(ChargeDensityData(:,1) / cfg.Ly);
    y_coor = mod(ChargeDensityData(:,1), cfg.Ly);
    ChargeDensity = (ChargeDensityData(:,2));
    SpinDensity = (SpinData(:,2));
end

% Calculate hole density
hole_density = 1 - ChargeDensity;

% Reshape hole density and spin data into matrices
hole_matrix = reshape(hole_density, [cfg.Ly, cfg.Lx]);
spin_matrix = reshape(SpinDensity, [cfg.Ly, cfg.Lx]);

% Create x and y coordinates for plotting (Meshgrid)
[x, y] = meshgrid(1:cfg.Lx, 1:cfg.Ly);
x = x(:);
y = y(:);

% Flatten hole and spin density for scatter and arrow plot
hole_flat = hole_matrix(:);
spin_flat = spin_matrix(:);
ave_hole = mean(hole_flat);
fprintf("Average hole density: %.4f\n", ave_hole);

% Normalize density to determine circle sizes
max_hole_val = max(hole_flat);
if max_hole_val == 0, max_hole_val = 1; end % Avoid div by zero
sizes = max_size * hole_flat / max_hole_val; 
sizes(sizes <= 0) = 1e-10; % Prevent zero/negative size errors

% ========== Plotting ==========
figure;
set(gcf, 'Color', 'white');
hold on;

% Draw only the outer boundary lines (Visual Box)
% Horizontal borders
plot([1, cfg.Lx], [1, 1], 'k', 'LineWidth', 2);
plot([1, cfg.Lx], [cfg.Ly, cfg.Ly], 'k', 'LineWidth', 2);
% Vertical borders
plot([1, 1], [1, cfg.Ly], 'k', 'LineWidth', 2);
plot([cfg.Lx, cfg.Lx], [1, cfg.Ly], 'k', 'LineWidth', 2);

% Define marker colors
marker_color = [190, 184, 220] / 220;

% Plot hole density as circles
scatter(x, y, sizes, 'filled', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', marker_color);

% Plot spin density as arrows
max_spin = max(abs(spin_flat));
if max_spin == 0, max_spin = 1; end

for i = 1:length(spin_flat)
    val = spin_flat(i);
    if val == 0, continue; end
    
    arrow_length = max_arrow_length * abs(val) / max_spin;
    arrow_shift = arrow_length / 2;  % Shift to center arrow in circle
    
    if val > 0
        % Up spin (Blue) - Shifted to visually center on point
        quiver(x(i), y(i) + arrow_shift, 0, -arrow_length, 0, 'Color', 'b', 'MaxHeadSize', 1, 'LineWidth', 1.5, 'AutoScale', 'off');
    elseif val < 0
        % Down spin (Red)
        quiver(x(i), y(i) - arrow_shift, 0, arrow_length, 0, 'Color', 'r', 'MaxHeadSize', 1, 'LineWidth', 1.5, 'AutoScale', 'off');
    end
end

% Set plot limits and axis properties
xlim([0.5, cfg.Lx + 3]);
ylim([0.5, cfg.Ly + 0.5]);
axis equal;
set(gca, 'YDir', 'reverse'); % Ensure origin is at the bottom-left (matrix coordinates)
axis off; 

% ========== Legend ==========
legend_x = cfg.Lx + 1;  % Position to the right of the plot
legend_y_start = cfg.Ly / 2 + 1;  % Start position
legend_spacing = 0.5;  % Vertical spacing

% Hole density legend
scatter(legend_x, legend_y_start, max_size * ave_hole / max_hole_val, 'filled', ...
    'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', marker_color);
text(legend_x + 0.3, legend_y_start, sprintf('Hole density = %.4f', ave_hole), 'FontSize', 14);

% Spin legend calculation
max_spin_val = max(abs(spin_flat));
if max_spin_val > 0
    demonstrate_spin = round(max_spin_val, digits - floor(log10(max_spin_val)) - 1);
    demo_length = max_arrow_length * demonstrate_spin / max_spin_val;
    
    % Spin Up example
    quiver(legend_x, legend_y_start - legend_spacing + demo_length/2, ...
        0, -demo_length, ...
        0, 'Color', 'b', 'MaxHeadSize', 1, 'LineWidth', 1.5, 'AutoScale', 'off');
    text(legend_x + 0.3, legend_y_start - legend_spacing, ...
        ['Spin Up = ', num2str(demonstrate_spin)], 'FontSize', 14);

    % Spin Down example
    quiver(legend_x, legend_y_start - 2 * legend_spacing - demo_length/2, ...
        0, demo_length, ...
        0, 'Color', 'r', 'MaxHeadSize', 1, 'LineWidth', 1.5, 'AutoScale', 'off');
    text(legend_x + 0.3, legend_y_start - 2 * legend_spacing, ...
        ['Spin Down = ', num2str(demonstrate_spin)], 'FontSize', 14);
end

hold off;
