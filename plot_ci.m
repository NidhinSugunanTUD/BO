%% compute mean, min, max
meanmetric_low = mean(all_metrics);
minmetric = min(all_metrics);
maxmetric = max(all_metrics);

x=1:51;
xconf = [x, x(end:-1:1)];
yconf = [maxmetric, minmetric(end:-1:1)];

figure(99)
p = fill(xconf, yconf, 10,"MarkerEdgeColor","#D95319", ...
    "MarkerFaceColor","#D95319",'FaceAlpha',0.2);

hold on;
plot(x, meanmetric_low, 'black')
hold off;

%%
meanmetric_low = mean(low_metrics);
meanmetric_med = mean(med_metrics);

std_dev_low = std(low_metrics, 0, 1);
std_dev_med = std(med_metrics, 0, 1);

curve_low1 = meanmetric_low + std_dev_low;
curve_low2 = meanmetric_low - std_dev_low;
curve_med1 = meanmetric_med + std_dev_med;
curve_med2 = meanmetric_med - std_dev_med;

x=1:51;
x2 = [x, fliplr(x)];
inBetween_low = [curve_low1, fliplr(curve_low2)];
inBetween_med = [curve_med1, fliplr(curve_med2)];

fill_color_low = [1, 0, 0];  % Red color
fill_color_med = [0, 0, 1];  % Blue color

fill(x2, inBetween_low, fill_color_low, 'FaceAlpha', 0.4, 'DisplayName', '');
hold on;
fill(x2, inBetween_med, fill_color_med, 'FaceAlpha', 0.1, 'DisplayName', '')
hold on;
plot(x, meanmetric_low, 'r', 'LineWidth', 2, 'DisplayName', 'Low noise');
hold on;
plot(x, meanmetric_med, 'b', 'LineWidth', 2, 'DisplayName', 'High noise');
hold on;
plot(x, zero_metrics, 'g', 'LineWidth', 2, 'DisplayName', 'No noise');
xlabel('Iteration');
ylabel('Objective function value');
legend('Location','best');
