close all; clear; clc
load ('A:\Matlab\MPC\BO_v_1\Algo noise\All metric\allmetric.mat')

zero_SDiter = mean(zero);
low_SDiter = std(low);
medium_SDiter = std(medium);
high_SDiter = std(high);

plot_allmetric = [low_SDiter; medium_SDiter; high_SDiter];
% plot_allmetric = [low_SDiter; medium_SDiter];

figure(501)
plot(zero_SDiter)
hold on
% bar(plot_allmetric')
bar(plot_allmetric','stacked')
legend('zero', 'low', 'medium', 'high')
% legend('zero', 'low', 'medium')
% legend('low', 'medium')
%%
close all; clear; clc
load ('A:\Matlab\MPC\BO_v_1\Algo noise\All metric\allmetric1.mat')
zero_SDiter = mean(zero);
low_SDiter = std(low);
medium_SDiter = std(medium);
high_SDiter = std(high1);

plot_allmetric = [low_SDiter; medium_SDiter; high_SDiter];
% plot_allmetric = [low_SDiter; medium_SDiter];

figure(501)
plot(zero_SDiter)
hold on
% bar(plot_allmetric')
bar(plot_allmetric','stacked')
legend('zero', 'low', 'medium', 'high')
% legend('zero', 'low', 'medium')
% legend('low', 'medium')
%%
close all; clear; clc
load ('A:\Matlab\MPC\BO_v_1\Algo noise\All metric\allmetric2.mat')
zero_SDiter = mean(zero);
low_SDiter = std(low);
medium_SDiter = std(medium);
high_SDiter = std(high2);

plot_allmetric = [low_SDiter; medium_SDiter; high_SDiter];
% plot_allmetric = [low_SDiter; medium_SDiter];

figure(501)
% plot(zero_SDiter)
% hold on
bar(plot_allmetric')
% bar(plot_allmetric','stacked')
legend('low', 'medium', 'high')
% legend('zero', 'low', 'medium')
% legend('low', 'medium')