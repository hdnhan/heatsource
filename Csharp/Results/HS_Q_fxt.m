clear
clc

nn = 64;
f = @(x, y, t) sin(pi * x).*y.*(1-y) * (t .^ 2 + 1);
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
legend('exact f(x_p, t)', 'f^\gamma_h with noise 1%', 'Location', 'south');
xlim([0, 1]); ylim ([0, 0.6]);
title(''); xlabel(''); ylabel(''); box on;
%xticks([0:0.1:1]);
saveas(gcf, 'HS_Q_fxt.png');