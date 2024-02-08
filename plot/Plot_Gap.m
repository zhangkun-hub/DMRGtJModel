Ly = [4, 6, 8];
gap = [0.21, 0.2, 0.5];

plot(Ly, gap, 'diamond'); hold on;


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2.5); % Set line width 1.5 pounds
xlabel('$L_y$','Interpreter','latex');
ylabel('SC gap','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

set(gca, 'XTick', [4,6,8]);
set(gca,'XTickLabel',{'4','6','8'});

set(gca, 'Ylim', [0, 0.6]);

set(gcf,'position',[1000,1000,450,350]);

figure_directory = './';
figure_name_eps = 'energy_gap.eps';
figure_path = fullfile(figure_directory, figure_name_eps);
saveas(gcf, figure_path, 'epsc');
disp(['energy gap figure saved as: ', figure_path]);