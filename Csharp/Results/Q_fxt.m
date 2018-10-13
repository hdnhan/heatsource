clear
clc

nn = 32;
f = @(x, y, t) (x .^ 3 + y .^ 3) * (t .^ 2 + 1);
time = 0:1/nn:1;

formatSpec = '%f';
fileID = fopen(strcat('Q_fxt_', num2str(nn), '_1%.txt'), 'r');
fh = fscanf(fileID, formatSpec);

fht = zeros(nn + 1, 1);
for nt = 0:nn
    for ny = 0:nn
        for nx = 0:nn
            if nx / nn == 0.5 && ny / nn == 0.5
                fht(nt + 1) = fh(nt * (nn + 1) ^ 2 + ny * (nn + 1) + nx + 1);
            end
        end
    end
end

figure
hold on
set(fplot(@(t)f(0.5, 0.5, t)), 'color', 'black')
plot(time, fht, 'black-o', 'MarkerSize', 2);
legend('exact f(0.5, 0.5, t)', 'f^\gamma_h with noise 1%', 'Location', 'south');
xlim([0, 1]); ylim auto;
title(''); xlabel(''); ylabel(''); box on;
%xticks([0:0.1:1]);