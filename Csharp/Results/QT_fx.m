clear
clc

nn = 32;
f = @(x, y) sin(pi*x).*y.*(1-y);
time = 0:1/nn:1;

formatSpec = '%f';
fileID = fopen(strcat('Q_fxt_', num2str(nn), '_1%.txt'), 'r');
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
xlim([0, 1]); ylim auto;
title(''); xlabel(''); ylabel(''); box on;

figure
hold on
set(fplot(@(y)f(0.5, y)), 'color', 'black')
plot(time, fhy, 'black-o', 'MarkerSize', 2);
legend('exact f(0.5, y)', 'f^\gamma_h with noise 1%', 'Location', 'south');
xlim([0, 1]); ylim auto;
title(''); xlabel(''); ylabel(''); box on;

%xticks([0:0.1:1]);