% Plot_SC_corr1D.m
% Plot 1D cut of Superconducting (SC) correlation function (log-log plot)
% Supports both OBC and Cylinder (PBC) geometries
% Replaces PlotSC_corr1D.m

clear;

% ====== Configuration Section ======
% Basic parameters from GetCommonConfig
cfg = GetCommonConfig();

% Visualization Parameters
legend_location = 'SouthWest';
font_size = 24;
line_width = 1.5;

% ===================================

figure;
lattice = SquareLattice(cfg.Ly, cfg.Lx);

[refer_x, ~] = lattice.indexToCoord(cfg.refer_site);

for i = 1:numel(cfg.D_values)
    Db = cfg.D_values(i);
    
    % Generate Filename Postfix
    StandardPostfix = GetSimulationParamsPostfix(cfg.Ly, cfg.Lx, cfg.t2, cfg.J, cfg.FillingType, cfg.FillingValue, Db, cfg.Geometry);
    FileNamePostfix = ['_fix_refer', num2str(cfg.refer_site), StandardPostfix];
    
    try
        % 1. Load and Plot Y-Y Bond Correlations ('scsyy')
        FileNamePrefix = '../data/scsyy';
        [sites_set, ~, sc_correlation] = LoadSCData(FileNamePrefix, FileNamePostfix, cfg.Ly, cfg.pairing_type);
        [target_x, ~] = lattice.indexToCoord(sites_set(:, 3));
        
        % Plot log-log
        loglog(target_x - refer_x, abs(sc_correlation), '-o', ...
            'DisplayName', ['D=', num2str(Db), ', yy'], ...
            'LineWidth', line_width);
        hold on;
        
        % 2. Load and Plot Y-X Bond Correlations ('scsyx')
        FileNamePrefix = '../data/scsyx';
        [sites_set, ~, sc_correlation] = LoadSCData(FileNamePrefix, FileNamePostfix, cfg.Ly, cfg.pairing_type);
        [target_x, ~] = lattice.indexToCoord(sites_set(:, 3));
        
        % Plot log-log
        loglog(target_x - refer_x, abs(sc_correlation), '-x', ...
            'DisplayName', ['D=', num2str(Db), ', yx'], ...
            'LineWidth', line_width);
            
    catch ME
        warning('Failed to load or plot data for D=%d: %s', Db, ME.message);
        continue;
    end
end

% Legend and Axis Styling
l = legend('Location', 'best');
set(l, 'Box', 'off', 'FontName', 'Arial', 'FontSize', 16, 'Location', legend_location);

xlabel('distance r', 'FontName', 'Arial');

if cfg.pairing_type == 's'
    ylabel('singlet-pair correlation', 'FontName', 'Arial');
elseif cfg.pairing_type == 't'
    ylabel('triplet-pair correlation', 'FontName', 'Arial');
end

set(gca, 'FontSize', font_size, 'LineWidth', line_width);
set(get(gca, 'Children'), 'LineWidth', line_width);
set(get(gca, 'XLabel'), 'FontSize', font_size);
set(get(gca, 'YLabel'), 'FontSize', font_size);
box on;
hold off;

