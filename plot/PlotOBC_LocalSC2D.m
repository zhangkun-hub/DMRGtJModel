% PlotOBC_LocalSC2D.m
clear;

% ====== Configuration ======
cfg = GetCommonConfig();

% Visualization Parameters
max_sc_bond_width = 20;
bond_width = 0.2;
bond_color = 'k';
positive_sc_color = [0     172	156]/256;
negative_sc_color = [237   159	84 ]/256;
reference_point_color = [152  115	232]/256;

% Loop over all bond dimensions in config
for i = 1:numel(cfg.D_values)
    Db = cfg.D_values(i);
    
    % Create a new figure for each Bond Dimension to avoid overplotting
    figure;
    set(gcf, 'Name', ['D = ' num2str(Db)]);
    
    lattice = SquareLattice(cfg.Ly, cfg.Lx);
    % sub_lattice = SquareLattice(cfg.Ly, floor(cfg.Lx/4)+3); % Unused or specific logic
    lattice.drawBonds(bond_color, bond_width);
    
    if(cfg.refer_site < 0) % Basic sanity check
        fprintf("Parameter refer site setting is suspicious: %d\n", cfg.refer_site);
    end
    
    % Create the file path
    % Note: This assumes filenames follow the _fix_referX pattern
    FileNamePostfix = ['_fix_refer', num2str(cfg.refer_site), ...
        GetSimulationParamsPostfix(cfg.Ly, cfg.Lx, cfg.t2, cfg.J, ...
        cfg.FillingType, cfg.FillingValue, Db, cfg.Geometry)];
        
    % Load Y-Y Data
    FileNamePrefixYY = '../data/scsyy';
    try
        [sites_set, ~, sc_correlation_yy] = LoadLocalSCData(FileNamePrefixYY, FileNamePostfix, cfg.Ly, cfg.pairing_type);
    catch ME
        warning('Skipping D=%d: Failed to load YY data: %s', Db, ME.message);
        close(gcf);
        continue;
    end
    
    refer_bond = sites_set(1, 1:2);
    target_bonds = sites_set(:, 3:4);
    
    % Scale width
    max_val = max(abs(sc_correlation_yy(:,1)));
    if max_val == 0, max_val = 1; end
    sc_bond_width_scale = max_sc_bond_width / max_val;
    
    % Plot y-y bond SC
    lattice.drawBondSC(target_bonds, sc_correlation_yy(:,1), sc_bond_width_scale, positive_sc_color, negative_sc_color);
    hold on;

    % Load Y-X Data
    FileNamePrefixYX = '../data/scsyx';
    try
        [sites_set, ~, sc_correlation_yx] = LoadLocalSCData(FileNamePrefixYX, FileNamePostfix, cfg.Ly, cfg.pairing_type);
        % Plot y-x bond SC
        lattice.drawBondSC(sites_set(:, 3:4), sc_correlation_yx(:,1), sc_bond_width_scale, positive_sc_color, negative_sc_color);
    catch
        warning('YX data missing for D=%d, skipping YX bonds.', Db);
    end

    % Legend
    lattice.drawOneBondData((cfg.Lx)*cfg.Ly+cfg.Lx/2, (cfg.Lx+1)*cfg.Ly+cfg.Lx/2, ...
        max_val, sc_bond_width_scale, ...
        positive_sc_color, negative_sc_color);
        
    text(cfg.Lx + 1, cfg.Ly/2, ...
        num2str(max_val), ...
        'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold');

    % Set plot limits and axis properties
    xlim([0.5, cfg.Lx + 5]);
    ylim([0.5, cfg.Ly + 0.5]);
    axis equal;
    set(gca, 'YDir', 'reverse'); % Ensure origin is at the bottom-left
    axis off; % Remove axis ticks and labels
    hold off;
    
    title(sprintf('Local SC (D=%d, Ref=%d)', Db, cfg.refer_site));
end
