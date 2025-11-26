% plot the charge density for cylinder systems.
% We assume the translational invariant of charge density in y direction, 
% and plot the charge density along x direction.
% The charge density is calculated by the formula:
% n(x) = 1/L_y \sum_{y} \langle n_{x, y} \rangle

clear;
Ly = 6; 
Lx = 48;
t2 = 0.2;
J = 0.33;
Hole = 48;
% Mu = 2;
D_values = [15000,20000,25000];
Geometry = 'Cylinder'; % Geometry type

D = D_values(1);
legendEntries = cell(1, numel(D_values));

if exist("Mu", "var")
    FillingType = 'Mu'; FillingValue = Mu;
elseif exist("Hole", "var")
    FillingType = 'Hole'; FillingValue = Hole;
end


ChargeDensityData = jsondecode(fileread(['../data/nf', GetSimulationParamsPostfix(Ly, Lx, t2, J, FillingType, FillingValue, D, Geometry)]));
position = zeros(1, size(ChargeDensityData,1));

if iscell(ChargeDensityData) % complex number DMRG
    for i=1:numel(position)
        position(i) = fix(ChargeDensityData{i}{1} / Ly);
    end
else % real number DMRG
    position = fix(ChargeDensityData(:,1)'/Ly);
end


ChargeDensity = zeros( numel(D_values), numel(position) );

for j = 1:numel(D_values)
    D=D_values(j);
    ChargeDensityData = jsondecode(fileread(['../data/nf', GetSimulationParamsPostfix(Ly, Lx, t2, J, FillingType, FillingValue, D, Geometry)]));
    if iscell(ChargeDensityData) % complex number DMRG
        for i=1:numel(position)
            ChargeDensity(j, i) = transpose(ChargeDensityData{i}{2}(1));
        end
    else % real number DMRG
        ChargeDensity(j,:) = ChargeDensityData(:,2)';
    end
    plot(position, ChargeDensity(j,:),'-x'); hold on;
    % print the info of average charge density
    fprintf('D=%d, <n>=%f\n', D, mean(ChargeDensity(j,:)));
    legendEntries{j} = ['$D=', num2str(D),'$'];
end
l = legend(legendEntries, 'Location', 'best', 'Interpreter', 'latex');
set(l, 'box', 'off');
set(l, 'fontsize', 24);

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('Charge Density','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
