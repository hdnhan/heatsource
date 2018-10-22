clear
clc

nn = 64;
f = @(t) 2 * t .* (t < 0.5) + 2 * (1 - t) .* (t >= 0.5);
time = 0:1/nn:1;

formatSpec = '%f';
fileID = fopen(strcat('Integration_ft_', num2str(nn), '_1%_1.txt'), 'r');
fh = fscanf(fileID, formatSpec);

figure
hold on
set(fplot(f), 'color', 'black');
plot(time, fh, 'black-o', 'MarkerSize', 2);
legend('exact f(t)', 'f^\gamma_h with noise 1%', 'Location', 'south');
xlim([0, 1]); ylim ([0, 1.2]);
title(''); xlabel(''); ylabel(''); box on;
%xticks([0:0.1:1]);
saveas(gcf, 'HS_Integration_ft1.png');



f = @(t) (t >= 0.25) .* (t <= 0.75);
fileID = fopen(strcat('Integration_ft_', num2str(nn), '_1%_2.txt'), 'r');
fh = fscanf(fileID, formatSpec);

figure
hold on
set(fplot(f), 'color', 'black');
plot(time, fh, 'black-o', 'MarkerSize', 2);
legend('exact f(t)', 'f^\gamma_h with noise 1%', 'Location', 'south');
xlim([0, 1]); ylim ([-0.2, 1.2]);
title(''); xlabel(''); ylabel(''); box on;
%xticks([0:0.1:1]);
saveas(gcf, 'HS_Integration_ft2.png');

