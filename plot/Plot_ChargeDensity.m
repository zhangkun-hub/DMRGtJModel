%

Lx=24; Ly=6;
t = 3;
t2 = 0;
J = 1;
J2 = 0;
Numhole = 6;

Dset = [10000];
D = Dset(1);
%nftJ8x24t2-0.6J20.04Hole32D15000.json
FileNamePostfix=['tJ',num2str(Ly),'x',num2str(Lx),'t2', num2str(t2),'J2', num2str(J2), 'Hole',num2str(Numhole),'D',num2str(D),'.json'];
ChargeDensityData = jsondecode(fileread(['../data/nf',FileNamePostfix]));
position = zeros(1, size(ChargeDensityData,1));
for i=1:numel(position)
    position(i) = fix((ChargeDensityData(i, 1))/Ly);
end

ChargeDensity = zeros( numel(Dset), numel(position) );

for j = 1:numel(Dset)
    D=Dset(j);
    FileNamePostfix=['tJ',num2str(Ly),'x',num2str(Lx),'t2', num2str(t2),'J2', num2str(J2), 'Hole',num2str(Numhole),'D',num2str(D),'.json'];
    ChargeDensityData = jsondecode(fileread(['../data/nf',FileNamePostfix]));
    ChargeDensity(j, :) = transpose(ChargeDensityData(:,2));
    plot(position, ChargeDensity(j,:),'-x'); hold on;
end

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('Charge Density','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
