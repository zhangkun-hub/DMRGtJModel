% Plot_Hopping_2D.m
% Visualizes the hopping (kinetic energy) terms on NN and NNN bonds.
% Checks for spin inversion symmetry (Up vs Down) and plots the sum.

clear;
close all;

% ====== Configuration Section ======
% Basic parameters are loaded from GetCommonConfig.m
cfg = GetCommonConfig();

% Override specific params if needed (or just use cfg directly)
% Local visualization params
bond_width_scale = 200; 
Theme = ModernPalettes.Get('AmberCyan');
pos_color = Theme.Positive;
neg_color = Theme.Negative;
check_tolerance = 1e-2;

% ===================================

% Generate Filename Postfix using config
Postfix = GetSimulationParamsPostfix(cfg.Ly, cfg.Lx, cfg.t2, cfg.J, cfg.FillingType, cfg.FillingValue, cfg.D_values(1), cfg.Geometry);

% ... (rest of logic uses cfg.Ly, cfg.Lx, etc.)
BondTypes = {
    'NN', 'hopping_nn_up', 'hopping_nn_down';
    'NNN', 'hopping_nnn_up', 'hopping_nnn_down'
};

figure;
lattice = SquareLattice(cfg.Ly, cfg.Lx);
lattice.drawBonds([0.8 0.8 0.8], 0.1); 
title_str = sprintf('Hopping Terms (Up+Down), %s', Postfix(1:end-5));
title(title_str, 'Interpreter', 'none');

% Iterate
for i = 1:size(BondTypes, 1)
    TypeLabel = BondTypes{i, 1};
    FilePrefixUp = BondTypes{i, 2};
    FilePrefixDown = BondTypes{i, 3};
    
    FileUp = ['../data/', FilePrefixUp, Postfix];
    FileDown = ['../data/', FilePrefixDown, Postfix];
    
    fprintf('Processing %s bonds...\n', TypeLabel);
    
    try
        DataUp = LoadHoppingData(FileUp);
        DataDown = LoadHoppingData(FileDown);
    catch ME
        warning('Skipping %s: %s', TypeLabel, ME.message);
        continue;
    end
    
    if numel(DataUp) ~= numel(DataDown)
        warning('Mismatch in data length for %s (Up: %d, Down: %d)', TypeLabel, numel(DataUp), numel(DataDown));
        continue;
    end
    
    MaxVal = 0;
    
    for k = 1:numel(DataUp)
        % Data Processing ... 
        EntryUp = DataUp{k};
        SitesUp = EntryUp{1}; 
        ValUp = EntryUp{2};
        
        EntryDown = DataDown{k};
        ValDown = EntryDown{2};
        
        % Symmetry Check
        Diff = abs(ValUp - ValDown);
        MeanVal = 0.5 * abs(ValUp + ValDown);
        if MeanVal > 1e-8 && (Diff / MeanVal > check_tolerance)
             warning('Spin symmetry violation at %s bond %d-%d: Up=%g, Down=%g, RelErr=%g', ...
                    TypeLabel, SitesUp(1), SitesUp(2), ValUp, ValDown, Diff/MeanVal);
        end
        
        TotalVal = ValUp + ValDown;
        lattice.drawOneBondData(SitesUp(1), SitesUp(2), TotalVal, bond_width_scale, pos_color, neg_color);
        MaxVal = max(MaxVal, abs(TotalVal));
    end
    fprintf('Finished %s bonds. Max magnitude: %g\n', TypeLabel, MaxVal);
end

% Adjustments
axis equal;
xlim([0.5, cfg.Lx + 1]);
ylim([0.5, cfg.Ly + 0.5]);
set(gca, 'YDir', 'reverse');
axis off;

% Helper
function Data = LoadHoppingData(Filename)
    if ~isfile(Filename)
        error('File not found: %s', Filename);
    end
    Raw = fileread(Filename);
    Data = jsondecode(Raw);
    if ~iscell(Data)
        error('Unexpected JSON format in %s', Filename);
    end
end

