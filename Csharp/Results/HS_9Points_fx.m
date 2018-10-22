clear
clc

figure
box on
xb = [0.1, 0.2, 0.4, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
yb = [0.4, 0.7, 0.1, 0.4, 0.9, 0.6, 0.2, 0.8, 0.5];
plot(xb, yb, 'black*')
xlim([0, 1]);
ylim([0, 1]);
grid on
ax = gca;
ax.GridLineStyle = '--';
ax.GridColor = 'black';
ax.LineWidth = 1;
saveas(gcf, 'Nine_Points.png')

nn = 64;
f = @(x, y) sin(pi * x).*y.*(1-y);
time = 0:1/nn:1;

formatSpec = '%f';
fileID = fopen(strcat('QT_fx_', num2str(nn), '_1%.txt'), 'r');
fh = fscanf(fileID, formatSpec);

fhx = zeros(nn + 1, 1);
fhy = zeros(nn + 1, 1);

for ny = 0:nn
    for nx = 0:nn
        if nx / nn == 0.5
            fhy(ny + 1) = fh(ny * (nn + 1) + nx + 1);
        end
         if ny / nn == 0.5
            fhx(nx + 1) = fh(ny * (nn + 1) + nx + 1);
        end
    end
end


figure
hold on
set(fplot(@(x)f(x, 0.5)), 'color', 'black')
plot(time, fhy, 'black-o', 'MarkerSize', 2);
legend('exact f(x, 0.5)', 'f^\gamma_h with noise 1%', 'Location', 'south');
xlim([0, 1]); ylim ([0, 0.3]);
title(''); xlabel(''); ylabel(''); box on;
saveas(gcf, 'HS_9Points_fx1.png');

figure
hold on
set(fplot(@(y)f(0.5, y)), 'color', 'black')
plot(time, fhy, 'black-o', 'MarkerSize', 2);
legend('exact f(0.5, y)', 'f^\gamma_h with noise 1%', 'Location', 'south');
xlim([0, 1]); ylim ([0, 0.3]);
title(''); xlabel(''); ylabel(''); box on;
saveas(gcf, 'HS_9Points_fx2.png');
%xticks([0:0.1:1]);