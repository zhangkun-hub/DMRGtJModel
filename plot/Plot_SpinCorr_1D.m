% Plot_SpinCorr_1D.m
clear;

% ====== Configuration ======
cfg = GetCommonConfig();
% Use the first D value if multiple are provided
Db = cfg.D_values(1); 

% ====== Load Data ======
FileNamePostfix = GetSimulationParamsPostfix(cfg.Ly, cfg.Lx, cfg.t2, cfg.J, cfg.FillingType, cfg.FillingValue, Db, cfg.Geometry);

try
    % Construct filenames dynamically based on reference site
    % e.g., sz0sz -> sz{refer_site}sz
    file_zz = ['../data/sz', num2str(cfg.refer_site), 'sz', FileNamePostfix];
    file_pm = ['../data/sp', num2str(cfg.refer_site), 'sm', FileNamePostfix];
    file_mp = ['../data/sm', num2str(cfg.refer_site), 'sp', FileNamePostfix];
    
    SpinDataZZ = jsondecode(fileread(file_zz));
    SpinDataPM = jsondecode(fileread(file_pm));
    SpinDataMP = jsondecode(fileread(file_mp));
catch ME
    error('Failed to load data (File: %s): %s', file_zz, ME.message);
end

% Extract x and y coordinates, charge density, and spin data
num_sites = size(SpinDataZZ, 1);
x_coor = zeros(1, num_sites);
y_coor = zeros(1, num_sites);
SpinDensity = zeros(1, num_sites);

% Parse Data (Handle both Cell Array and Matrix formats)
if iscell(SpinDataZZ) 
    % Real number DMRG or specific JSON format (Cell Array)
    for i = 1:num_sites
        % SpinDataZZ{i}{1} is typically [ref_site, target_site]
        site_idx = SpinDataZZ{i}{1}(2);
        x_coor(i) = fix(site_idx / cfg.Ly) + 1;
        y_coor(i) = mod(site_idx, cfg.Ly) + 1;
        
        % Values are in the second element of the cell
        zz = SpinDataZZ{i}{2};
        pm = SpinDataPM{i}{2};
        mp = SpinDataMP{i}{2};
        
        SpinDensity(i) = zz + 0.5 * (pm + mp);
    end
else 
    % Complex number DMRG or Matrix format
    % Assumes columns: [ref_site, val_real, target_site, val_imag]
    for i = 1:num_sites
        site_idx = SpinDataZZ(i, 3);
        x_coor(i) = fix(site_idx / cfg.Ly) + 1;
        y_coor(i) = mod(site_idx, cfg.Ly) + 1;
        
        % Calculate Spin Correlation: <S(0)S(r)> = <SzSz> + 0.5(<S+S-> + <S-S+>)
        szsz = SpinDataZZ(i, 2) + 1i * SpinDataZZ(i, 4);
        spsm = SpinDataPM(i, 2) + 1i * SpinDataPM(i, 4);
        smsp = SpinDataMP(i, 2) + 1i * SpinDataMP(i, 4);
        
        SpinDensity(i) = szsz + 0.5 * (spsm + smsp);
    end
end

% Plotting
figure;
% Example logic: Plotting log-log of first 7 points (or fewer if data is small)
range = 1:min(7, numel(y_coor)); 
loglog(y_coor(range)-1, abs(SpinDensity(range)), 'o');
hold on;

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2);
xlabel('$\Delta x$','Interpreter','latex');
ylabel('$\langle \mathbf S(0) \cdot \mathbf S(r)\rangle$','Interpreter','latex')
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

hold off;
