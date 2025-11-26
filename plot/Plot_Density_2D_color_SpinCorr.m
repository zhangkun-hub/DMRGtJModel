% Plot_Density_2D_color_SpinCorr.m
% Combined plot of Charge Density (colormap) and Spin Correlations (arrows)
% Refactored to use centralized configuration

clear;
close all;

% ====== Configuration ======
cfg = GetCommonConfig();

% Visualization Parameters
max_arrow_length = 0.6;
digits = 2;     % Significant figures for legend

% ===== IMPROVED COLOR SCHEME =====
% High-contrast colors for better visibility
positive_spin_color = [0.95, 0.45, 0.05];  % Bright gold (up arrows)
negative_spin_color = [0.45, 0.0, 0.75];    % Deep purple (down arrows)

% Use the first D value if multiple are provided
Db = cfg.D_values(1);

% Generate Filename Postfix
FileNamePostfix = GetSimulationParamsPostfix(cfg.Ly, cfg.Lx, cfg.t2, cfg.J, ...
    cfg.FillingType, cfg.FillingValue, Db, cfg.Geometry);

% ========== Load Hole Density Data ==========
try
    charge_density_file = ['../data/nf', FileNamePostfix];
    fprintf("Loading charge density data from %s\n", charge_density_file);
    ChargeDensityData = jsondecode(fileread(charge_density_file));
catch ME
    error('Failed to load charge density data: %s', ME.message);
end

% Extract coordinates and densities
if iscell(ChargeDensityData)
    num_points = numel(ChargeDensityData);
    x_coor = zeros(1, num_points);
    y_coor = zeros(1, num_points);
    ChargeDensity = zeros(1, num_points);
    for i = 1:num_points
        site_idx = ChargeDensityData{i}{1};
        x_coor(i) = fix(site_idx / cfg.Ly) + 1;
        y_coor(i) = mod(site_idx, cfg.Ly) + 1;
        ChargeDensity(i) = ChargeDensityData{i}{2}(1);
    end
else
    site_indices = ChargeDensityData(:,1);
    x_coor = fix(site_indices / cfg.Ly) + 1;
    y_coor = mod(site_indices, cfg.Ly) + 1;
    ChargeDensity = ChargeDensityData(:,2);
end

% Calculate hole density
hole_density = 1 - ChargeDensity;
hole_matrix = reshape(hole_density, [cfg.Ly, cfg.Lx]);
ave_hole_density = mean(hole_density);
fprintf("Average hole density: %.4f\n", ave_hole_density);

% ========== Load Spin Correlation Data ==========
try
    % Construct filenames dynamically based on reference site
    file_zz = ['../data/sz', num2str(cfg.refer_site), 'sz', FileNamePostfix];
    file_pm = ['../data/sp', num2str(cfg.refer_site), 'sm', FileNamePostfix];
    file_mp = ['../data/sm', num2str(cfg.refer_site), 'sp', FileNamePostfix];
    
    SpinDataZZ = jsondecode(fileread(file_zz));
    SpinDataPM = jsondecode(fileread(file_pm));
    SpinDataMP = jsondecode(fileread(file_mp));
catch ME
    error('Failed to load spin correlation data: %s', ME.message);
end

% Process spin correlation data
num_spin_points = size(SpinDataZZ, 1);
x_plot = zeros(num_spin_points, 1);
y_plot = zeros(num_spin_points, 1);
SpinCorr = zeros(num_spin_points, 1);

if iscell(SpinDataZZ)  % Real DMRG data
    for i = 1:num_spin_points
        site_idx = SpinDataZZ{i}{1}(2);
        x_plot(i) = fix(site_idx / cfg.Ly) + 1;
        y_plot(i) = mod(site_idx, cfg.Ly) + 1;
        
        zz_val = SpinDataZZ{i}{2}(1);
        pm_val = SpinDataPM{i}{2}(1);
        mp_val = SpinDataMP{i}{2}(1);
        
        SpinCorr(i) = zz_val + 0.5*(pm_val + mp_val);
    end
else  % Complex DMRG data
    for i = 1:num_spin_points
        site_idx = SpinDataZZ(i,3);
        x_plot(i) = fix(site_idx / cfg.Ly) + 1;
        y_plot(i) = mod(site_idx, cfg.Ly) + 1;
        
        zz_val = SpinDataZZ(i,2);
        pm_val = SpinDataPM(i,2);
        mp_val = SpinDataMP(i,2);
        
        SpinCorr(i) = zz_val + 0.5*(pm_val + mp_val);
    end
end

% ========== Plotting ==========
figure;
set(gcf, 'Position', [100, 100, 800, 500], 'Color', 'white');
set(gcf, 'CloseRequestFcn', @(~,~) closereq);

% Plot hole density colormap
imagesc(1:cfg.Lx, 1:cfg.Ly, hole_matrix);
set(gca, 'YDir', 'reverse');
try
    colormap(slanCM(19));  
catch
    colormap('parula');
end
cbar = colorbar;
ylabel(cbar, 'Hole Density', 'FontSize', 12);
hold on;

% ===== IMPROVED ARROW PLOTTING =====
% Pre-calculate arrow properties
max_spin = max(abs(SpinCorr));
if max_spin > 0
    arrow_scale = max_arrow_length / max_spin;
else
    arrow_scale = 0;
end

% Plot spin correlation arrows with enhanced appearance
for i = 1:length(SpinCorr)
    if SpinCorr(i) == 0, continue; end
    
    % Calculate arrow properties
    arrow_length = abs(SpinCorr(i)) * arrow_scale;
    arrow_head_size = 0.5 + 1.5 * (abs(SpinCorr(i)))/max_spin;  % Dynamic head size
    arrow_shift = arrow_length / 2;  % Shift to center arrow in circle
    
    if SpinCorr(i) > 0  % Up spin
        quiver(x_plot(i), y_plot(i) + arrow_shift, 0, -arrow_length, ...
            0, 'Color', positive_spin_color, 'LineWidth', 2.2, ...
            'MaxHeadSize', arrow_head_size, 'AutoScale', 'off', ...
            'AlignVertexCenters', 'on');
    else  % Down spin
        quiver(x_plot(i), y_plot(i) - arrow_shift, 0, arrow_length, ...
            0, 'Color', negative_spin_color, 'LineWidth', 2.2, ...
            'MaxHeadSize', arrow_head_size, 'AutoScale', 'off', ...
            'AlignVertexCenters', 'on');
    end
end

% Add reference site marker (more prominent)
ref_x = fix(cfg.refer_site / cfg.Ly) + 1;
ref_y = mod(cfg.refer_site, cfg.Ly) + 1;
plot(ref_x, ref_y, 'o', 'MarkerSize', 20, 'LineWidth', 2.5, ...
    'MarkerEdgeColor', 'k');

% Create legend
max_spin = max(abs(SpinCorr));
if max_spin > 0 && isfinite(max_spin)
    demo_value = round(max_spin, digits - floor(log10(max_spin)) - 1);
    demo_arrow_length = demo_value * arrow_scale;

    % Positive example
    quiver(cfg.Lx+0.8, cfg.Ly-1, 0, -demo_arrow_length, 0, 'Color', positive_spin_color, ...
        'LineWidth', 2.2, 'MaxHeadSize', 1.2, 'AutoScale','off', ...
        'AlignVertexCenters', 'on');
    text(cfg.Lx+1.5, cfg.Ly-1, sprintf('Corr = %.2f', demo_value), ...
        'FontSize', 11, 'VerticalAlignment','middle', 'FontWeight','bold');

    % Negative example
    quiver(cfg.Lx+0.8, cfg.Ly-2, 0, demo_arrow_length, 0, 'Color', negative_spin_color, ...
        'LineWidth', 2.2, 'MaxHeadSize', 1.2, 'AutoScale','off', ...
        'AlignVertexCenters', 'on');
    text(cfg.Lx+1.5, cfg.Ly-2, sprintf('Corr = -%.2f', demo_value), ...
        'FontSize', 11, 'VerticalAlignment','middle', 'FontWeight','bold');
end


% Final adjustments
xlim([0.5, cfg.Lx+3]);
ylim([0.5, cfg.Ly+0.5]);
axis equal tight;
set(gca, 'XTick', [], 'YTick', [], 'Box', 'off');
axis off;
drawnow;

% Add title with parameters (Optional, commented out in original)
% title_str = sprintf('Ly=%d, Lx=%d, J=%.2f, Filling=%.2f, D=%d', cfg.Ly, cfg.Lx, cfg.J, cfg.FillingValue, Db);
% title(title_str, 'FontSize', 12, 'Interpreter', 'none');

hold off;
