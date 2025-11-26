% Plot_SpinCorr_2D.m
% Plot 2D Spin Correlation functions (arrows)
% Supports both OBC and Cylinder (PBC) geometries
% Replaces PlotOBC_SpinCorr.m

clear;
close all;

% ====== Configuration Section ======
% Basic parameters from GetCommonConfig
cfg = GetCommonConfig();

% Visualization Parameters
max_arrow_length = 1;   % Max arrow length for visual scaling
use_cubic_root = false;  % Apply cubic root scaling to flatten data (visual trick)

% Colors
positive_spin_color = 'b';
negative_spin_color = 'r';
marker_colors = {[190, 184, 220] / 220}; % Light purple

% ===================================

% Generate Filename (Using first D_value from config if list provided)
Db = cfg.D_values(1);
FileNamePostfix = GetSimulationParamsPostfix(cfg.Ly, cfg.Lx, cfg.t2, cfg.J, cfg.FillingType, cfg.FillingValue, Db, cfg.Geometry);

try
    % Load Data
    SpinDataZZ = jsondecode(fileread(['../data/sz', num2str(cfg.refer_site), 'sz', FileNamePostfix]));
    SpinDataPM = jsondecode(fileread(['../data/sp', num2str(cfg.refer_site), 'sm', FileNamePostfix]));
    SpinDataMP = jsondecode(fileread(['../data/sm', num2str(cfg.refer_site), 'sp', FileNamePostfix]));
catch ME
    error('Failed to load spin data: %s', ME.message);
end

% Extract Data
% Determine if data is Real or Complex format based on JSON structure
is_complex = ~iscell(SpinDataZZ) && size(SpinDataZZ, 2) >= 4; 

num_sites = size(SpinDataZZ, 1);
x_coor = zeros(1, num_sites);
y_coor = zeros(1, num_sites);
SpinCorr = zeros(1, num_sites);

if iscell(SpinDataZZ) % Real number DMRG
    for i = 1:num_sites
        site_idx = SpinDataZZ{i}{1}(2);
        x_coor(i) = fix(site_idx / cfg.Ly) + 1;
        y_coor(i) = mod(site_idx, cfg.Ly) + 1;
        
        zz = SpinDataZZ{i}{2};
        pm = SpinDataPM{i}{2};
        mp = SpinDataMP{i}{2};
        
        SpinCorr(i) = zz + 0.5 * (pm + mp);
    end
else % Complex number DMRG (assumed structure)
    for i = 1:num_sites
        site_idx = SpinDataZZ(i, 3);
        x_coor(i) = fix(site_idx / cfg.Ly) + 1;
        y_coor(i) = mod(site_idx, cfg.Ly) + 1;
        
        zz = SpinDataZZ(i, 2) + 1i * SpinDataZZ(i, 4);
        pm = SpinDataPM(i, 2) + 1i * SpinDataPM(i, 4);
        mp = SpinDataMP(i, 2) + 1i * SpinDataMP(i, 4);
        
        SpinCorr(i) = real(zz + 0.5 * (pm + mp)); % Take real part for plotting
    end
end

% Check reference point consistency
% refer_site_check = SpinDataZZ(1, 1); 

% Visualization Processing
spin_flat = SpinCorr;
if use_cubic_root
    spin_flat = sign(SpinCorr) .* abs(SpinCorr).^(1/3);
    fprintf("Note: Spin magnitudes are cubic-root scaled for visualization!\n");
end

% Plotting
figure;
hold on;

% Plot Arrows
max_val = max(abs(spin_flat));
if max_val == 0, max_val = 1; end % Avoid div by zero

for i = 1:length(spin_flat)
    val = spin_flat(i);
    if val == 0, continue; end
    
    arrow_length = max_arrow_length * abs(val) / max_val;
    arrow_shift = arrow_length / 2;
    
    if val > 0
        % Up spin (Blue) - shifted up visually
        quiver(x_coor(i), y_coor(i) + arrow_shift, 0, -arrow_length, 0, ...
            'Color', positive_spin_color, 'MaxHeadSize', 1, 'LineWidth', 1.5);
    else
        % Down spin (Red) - shifted down visually
        quiver(x_coor(i), y_coor(i) - arrow_shift, 0, arrow_length, 0, ...
            'Color', negative_spin_color, 'MaxHeadSize', 1, 'LineWidth', 1.5);
    end
end

% Axis Settings
xlim([0.5, cfg.Lx + 3]);
ylim([0.5, cfg.Ly + 0.5]);
axis equal;
set(gca, 'YDir', 'reverse');
axis off;

% Legend
legend_x = cfg.Lx + 1;
legend_y_start = cfg.Ly / 2 + 1;
legend_spacing = 0.5;

digits = 1;
max_spin_real = max(abs(SpinCorr)); % Use real values for legend text
if max_spin_real > 0
    demo_spin = round(max_spin_real, digits - floor(log10(abs(max_spin_real))) - 1);
    
    % Calculate visual length for the demo value
    % Must apply same scaling logic!
    if use_cubic_root
        demo_val_scaled = abs(demo_spin)^(1/3);
    else
        demo_val_scaled = abs(demo_spin);
    end
    
    demo_length = max_arrow_length * demo_val_scaled / max_val;
    
    % Up Legend
    quiver(legend_x, legend_y_start - legend_spacing + demo_length/2, ...
        0, -demo_length, 0, 'Color', positive_spin_color, 'MaxHeadSize', 1, 'LineWidth', 1.5);
    text(legend_x + 0.3, legend_y_start - legend_spacing, ...
        ['Spin Up = ', num2str(demo_spin)], 'FontSize', 18);
        
    % Down Legend
    quiver(legend_x, legend_y_start - 2*legend_spacing - demo_length/2, ...
        0, demo_length, 0, 'Color', negative_spin_color, 'MaxHeadSize', 1, 'LineWidth', 1.5);
    text(legend_x + 0.3, legend_y_start - 2*legend_spacing, ...
        ['Spin Down = ', num2str(demo_spin)], 'FontSize', 18);
end

hold off;

