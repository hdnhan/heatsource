clear
clc

f = @(x) sin(pi * x);
x = [0:1/50:1];
figure
hold on
set(ezplot(f), 'color', 'black')
plot(x, f(x)+(0.5-rand(1, 51))/5, 'black-o', 'MarkerSize', 2);
legend('exact', 'approximate', 'Location','south');
xlim([0, 1]); ylim auto;
title(''); xlabel(''); ylabel(''); box on;
%xticks([0:0.1:1]);

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
