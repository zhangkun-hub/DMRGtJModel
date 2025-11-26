% Plot_SC_corr2D.m
% Unified script for plotting 2D Superconducting (SC) correlation functions
% Supports both OBC and Cylinder (PBC) geometries
% Replaces PlotOBC_SC_corr2D.m and PlotPBC_SC_corr2D.m

clear;

% ====== Configuration Section ======
% Basic parameters from GetCommonConfig
cfg = GetCommonConfig();

% Visualization Parameters
bond_width = 0.2;
bond_color = 'k';

% Choose a Theme from ModernPalettes
Theme = ModernPalettes.Get('EmeraldViolet'); 
positive_sc_color = Theme.Positive;
negative_sc_color = Theme.Negative;

reference_point_color = [152  115   232]/256; % Purple
max_sc_bond_width = 20; % Max width for strongest bond
draw_sublattice = false; % If true, only draws a portion of bonds (OBC style)
use_sqrt_scale = false; % If true, plot sqrt(abs(correlation)) to enhance visibility of weak signals

% ===================================

figure;
lattice = SquareLattice(cfg.Ly, cfg.Lx);

% Draw background lattice bonds
if draw_sublattice
    sub_lattice = SquareLattice(cfg.Ly, floor(cfg.Lx/4)+3);
    sub_lattice.drawBonds(bond_color, bond_width);
else
    lattice.drawBonds(bond_color, bond_width);
end

% Mark reference site
lattice.drawOneBondData(cfg.refer_site, cfg.refer_site+1, 3, 1, 'k', 'k');

for i = 1:numel(cfg.D_values)
    Db = cfg.D_values(i);
    
    % Generate Filename Postfix using standard helper
    StandardPostfix = GetSimulationParamsPostfix(cfg.Ly, cfg.Lx, cfg.t2, cfg.J, cfg.FillingType, cfg.FillingValue, Db, cfg.Geometry);
    FileNamePostfix = ['_fix_refer', num2str(cfg.refer_site), StandardPostfix];
    
    % 1. Load and Plot Y-Y Bond Correlations ('scsyy')
    FileNamePrefix = '../data/scsyy';
    try
        [sites_set, ~, sc_correlation] = LoadSCData(FileNamePrefix, FileNamePostfix, cfg.Ly, cfg.pairing_type);
        
        target_bonds = sites_set(:, 3:4);
        
        % Apply Sqrt Scaling if requested
        if use_sqrt_scale
            sc_correlation = sign(sc_correlation) .* sqrt(abs(sc_correlation));
        end

        % Determine scaling based on max value
        sc_bond_width_scale = max_sc_bond_width / max(abs(sc_correlation));
        
        lattice.drawBondSC(target_bonds, sc_correlation, sc_bond_width_scale, positive_sc_color, negative_sc_color);
        hold on;
        
        % 2. Load and Plot Y-X Bond Correlations ('scsyx')
        FileNamePrefix = '../data/scsyx';
        [sites_set, ~, sc_correlation_yx] = LoadSCData(FileNamePrefix, FileNamePostfix, cfg.Ly, cfg.pairing_type);
        target_bonds_yx = sites_set(:, 3:4);
        
        if use_sqrt_scale
            sc_correlation_yx = sign(sc_correlation_yx) .* sqrt(abs(sc_correlation_yx));
        end
        
        lattice.drawBondSC(target_bonds_yx, sc_correlation_yx, sc_bond_width_scale, positive_sc_color, negative_sc_color);
        
    catch ME
        warning('Failed to load or plot data for D=%d: %s', Db, ME.message);
        continue;
    end
end

% Legend Drawing
legend_x1 = (cfg.Lx)*cfg.Ly + cfg.Lx/2;
legend_x2 = (cfg.Lx+1)*cfg.Ly + cfg.Lx/2;

if exist('sc_correlation', 'var')
    max_corr = max(abs(sc_correlation));
    lattice.drawOneBondData(legend_x1, legend_x2, max_corr, sc_bond_width_scale, positive_sc_color, negative_sc_color);
    
    text(cfg.Lx + 1, cfg.Ly/2, ...
        num2str(max_corr, '%.2e'), ...
        'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold');
end

xlim([0.5, cfg.Lx + 5]);
ylim([0.5, cfg.Ly + 0.5]);
axis equal;
set(gca, 'YDir', 'reverse'); 
axis off;
hold off;

