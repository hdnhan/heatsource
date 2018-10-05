clear
clc
f = @(x, y, z) sin(pi * x) * (2 * y * ( y < 0.5) + 2 * (1 - y) * (y >= 0.5)) * ((z >= 0.25) * (z <= 0.75));

formatSpec = '%f';
% infor
fileID = fopen('2D_Q_fh.sol','r');
fh1 = fscanf(fileID, formatSpec);
fileID = fopen('2D_Integration_fh.sol','r');
fh2 = fscanf(fileID, formatSpec);
fileID = fopen('2D_Points_fh.sol','r');
fh3 = fscanf(fileID, formatSpec);

nn = 40;
k = 25;
for i=1:nn+1
    fh1x(i) = fh1(k*41*41+k*41+i);
    fh2x(i) = fh2(k*41*41+k*41+i);
    fh3x(i) = fh3(k*41*41+k*41+i);
end

for i=0:nn
    fh1y(i+1) = fh1(k*41*41+i*41+k);
    fh2y(i+1) = fh2(k*41*41+i*41+k);
    fh3y(i+1) = fh3(k*41*41+i*41+k);
end

for i=0:nn
    fh1t(i+1) = fh1(i*41*41+k*41+k);
    fh2t(i+1) = fh2(i*41*41+k*41+k);
    fh3t(i+1) = fh3(i*41*41+k*41+k);
end

figure
hold on
vari = [0:1/40:1];
plot(vari, fh1x, 'black->', 'MarkerSize', 3);
plot(vari, fh2x, 'red-o', 'MarkerSize', 3);
plot(vari, fh3x, 'green-+', 'MarkerSize', 3);
ezplot(@(x)f(x, 0.6, 0.6));
legend('fh Q', 'fh Integration', 'fh Point at (0.5, 0.5)', 'Exact f(x, 0.6, 0.6)', 'Location','south');
xlim ([0, 1])
ylim ([-0.1, 0.85])
title(''); xlabel(''); ylabel(''); box on;


figure
hold on
plot(vari, fh1y, 'black->', 'MarkerSize', 3);
plot(vari, fh2y, 'red-o', 'MarkerSize', 3);
plot(vari, fh3y, 'green-+', 'MarkerSize', 3);
ezplot(@(y)f(0.6, y, 0.6));
legend('fh Q', 'fh Integration', 'fh Point at (0.5, 0.5)', 'Exact f(0.6, y, 0.6)', 'Location','south');
xlim ([0, 1])
ylim ([-0.1, 1.05])
title(''); xlabel(''); ylabel(''); box on;


figure
hold on
plot(vari, fh1t, 'black->', 'MarkerSize', 3);
plot(vari, fh2t, 'red-o', 'MarkerSize', 3);
plot(vari, fh3t, 'green-+', 'MarkerSize', 3);
ezplot(@(t)f(0.6, 0.6, t));
legend('fh Q', 'fh Integration', 'fh Point at (0.5, 0.5)', 'Exact f(0.6, 0.6, t)', 'Location','south');
xlim ([0, 1])
ylim ([-0.1, 0.8])
title(''); xlabel(''); ylabel(''); box on;