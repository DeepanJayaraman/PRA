function A = Area_shade_plot(x,Act,E,cmom10,cmom50,cmom90,Clr,Legend)
% x - population 
% Act - Actual failures from weibul and log-weibull disributions
% E - Temperature 
% Cmom10 - 10th or 25th percentile vaules
% Cmom50 - 50th percentile values
% Cmom90 - 90th or 75th percentile vaules
% clr - specify C- or L-moment
figure;
switch Clr
    case 'C-moment'
        colour = 'k';
    case 'L-moment'
        colour = 'b';
end
facealpha = 0.05;

curve1 = cmom10;
curve2 = cmom90;

plot(x,Act,'r', 'LineWidth', 1.5); hold on; % Population

plot(E, curve1,colour , 'LineWidth', 1.5,'linestyle', ':'); % 10th or 25th percentile

plot(E,cmom50, colour,'linestyle', '--' ); % Median
hold on;

plot(E, curve2, colour, 'LineWidth', 1.5,'linestyle', '-.'); % 90th or 75th percentile

AF = [E;curve1']';
AA = flip(AF,1);

patch([E AA(:,1)'], [curve2; AA(:,2)], colour, 'facealpha',facealpha,'edgecolor','none')

% x2 = [E,E];
% inBetween = [fliplr(curve1); curve2];
% fill(x2, inBetween, colour);

xlim([42 70]);
ylim([-2 10]);
% set(gca,'yscale','log');

xlabel('Field joint temperature, {^\circF}')
ylabel('\it{E_{N_f}}')

legend(Legend)

set(gca,'fontname','times','fontsize',12,'fontweight','bold')

