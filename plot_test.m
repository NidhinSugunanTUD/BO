close; clear; clc;
% Define the quadratic equation
a_values = [0.1, 0.5, 1];  % Different values of 'a'
x1 = linspace(-10, 10, 1000);  % Generate x1 values
x2 = linspace(-10, 10, 1000);  % Generate x2 values

% Create a meshgrid of x1 and x2
[X1, X2] = meshgrid(x1, x2);

% Initialize a figure for plotting
figure;

% Plot the difference between x1 and x2 (error) with y as the x-axis
for a = a_values
    % Calculate the quadratic equation
    y = a * (X1 - X2).^2;

    % Calculate the error (difference between x1 and x2)
    error = X1 - X2;

    % Plot the error against y
    subplot(1, numel(a_values), find(a_values == a));
    plot(error(:), y(:));
    title(['a = ' num2str(a)]);
    xlabel('error');
    ylabel('y');
    
    % % Set y-axis limits
    % ylim([200, 0]);
    

end

% Adjust subplot layout
sgtitle('Difference between x1 and x2 vs. y');


%% standard deviation as a shaded area 
% https://nl.mathworks.com/matlabcentral/answers/494515-plot-standard-deviation-as-a-shaded-area
% color: https://nl.mathworks.com/help/matlab/ref/fill.html
y = rand(1,10); % your mean vector;
x = 1:numel(y);
std_dev = 0.5;
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 75,"MarkerEdgeColor",[0 0.4470 0.7410], ...
    "MarkerFaceColor",[0 0.4470 0.7410],'FaceAlpha',0.2);
hold on;
plot(x, y, 'b', 'LineWidth', 2);

%% Heatmap
% Sample data
data = sim.metric; % Replace this with your actual data

% Define x and y axis limits
x_limits = [0 50];
y_limits = [0 50];

% Scale the data to fit within the specified axis limits
scaled_data = data * (x_limits(2) - x_limits(1));

% Create a heatmap using imagesc
imagesc(linspace(x_limits(1), x_limits(2), size(data, 2)), ...
    linspace(y_limits(1), y_limits(2), size(data, 1)), scaled_data);

% Customize colormap
colormap(flipud(hot)); % 'hot' colormap goes from black to red, we flip it to start from red

% Add colorbar and adjust labels
colorbar('southoutside');

% Set axis limits
xlim(x_limits);
ylim(y_limits);

% Add labels and title
xlabel('X-axis label');
ylabel('Y-axis label');
title('Heatmap Title');

%% Objective Function Evolution plot
% Objective function values 
objective_value = sim.metric; 

% Number of iterations or evaluations
num_iterations = numel(objective_value);

% Plotting the Objective Function Evolution
figure;
plot(1:num_iterations, objective_value);
xlabel('Iterations');
ylabel('Objective Function Value');
title('Objective Function Evolution');
grid on;

%% Scatter Plot of Hyperparameter Combinations 
% Hyperparameter combinations and corresponding objective function values

Q_weights = [sim.W(1, :)',sim.W(2, :)'];
objective_values = sim.metric; 
BestQ1 = opt_weights(1);
BestQ2 = opt_weights(2);

% Plotting the Scatter Plot of Hyperparameter Combinations
figure;
scatter(Q_weights(:, 1), Q_weights(:, 2), 50, objective_values, 'filled');
colorbar; % Add colorbar to show objective function values
hold on;
scatter(BestQ1, BestQ2, 120, minBEST, 'rs', 'filled')

xlabel('Q_1');
ylabel('Q_4');
title('Scatter Plot of Hyperparameter Combinations');

%% R weights
% Hyperparameter combinations and corresponding objective function values

R_weights = [sim.W(1, :)',sim.W(3, :)'];
objective_values = sim.metric; 
BestR1 = opt_weights(1);
BestR2 = opt_weights(3);

% Plotting the Scatter Plot of Hyperparameter Combinations
figure;
scatter(R_weights(:, 1), R_weights(:, 2), 50, objective_values, 'filled');
colorbar; % Add colorbar to show objective function values
hold on;
scatter(BestR1, BestR2, 120, minBEST, 'rs', 'filled')

xlabel('R_1');
ylabel('R_2');
title('Scatter Plot of Hyperparameter Combinations');

%% Contour Plot of Objective Function
% Define the range of Q1 and Q2 values
Q1_values = linspace(min(sim.W(1, :)), max(sim.W(1, :)), 51); 
Q2_values = linspace(min(sim.W(2, :)), max(sim.W(2, :)), 51); 

% Create a grid of Q1 and Q2 values
[Q1_grid, Q2_grid] = meshgrid(Q1_values, Q2_values);

objective_values_grid = zeros(size(Q1_grid));
for i = 1:length(Q1_grid)
    objective_values_grid(i,:) = sim.metric(i);
end

% Create a contour plot
figure;
contourf(Q1_grid, Q2_grid, objective_values_grid);
colorbar; % Add colorbar to show objective function values
xlabel('Q1');
ylabel('Q2');
title('Contour Plot of Objective Function');

%% Sensitivity analysis

% Initialize matrix to store performance metric values
performance_metric_values = zeros(length(Q1_values), length(Q2_values));

% Evaluate performance metric for each combination of Q1 and Q2
for i = 1:length(Q1_values)
    for j = 1:length(Q2_values)
        Q_weights1 = [Q1_values(i), Q2_values(j)]; % Combine Q1 and Q2 into weights vector
        performance_metric_values(i, j) = sqrt(mean(error.^2));
    end
end

% Plot Sensitivity Analysis
figure;
contourf(Q1_values, Q2_values, performance_metric_values);
colorbar; % Add colorbar to show performance metric values
xlabel('Q1');
ylabel('Q2');
title('Sensitivity Analysis Plot of Performance Metric');

%% Heatmap

% % Initialize matrix to store performance metric values
% performance_metric_values = zeros(length(Q1_values), length(Q2_values));
% 
% % Evaluate performance metric for each combination of Q1 and Q2
% performance_metric_values = diag(sort_metric);
% 
% % Create Heatmap of Hyperparameter Performance
% figure;
% heatmap(Q1_values, Q2_values, performance_metric_values);
% xlabel('Q1');
% ylabel('Q2');
% title('Heatmap of Hyperparameter Performance');

%% 
best1 = zeros(1,length(sim1.metric));
best2 = zeros(1,length(sim2.metric));
best3 = zeros(1,length(sim3.metric));

for i = 1:length(sim1.metric)
    best1(i) = min(sim1.metric(1:i));
end

for j = 1:length(sim2.metric)
    best2(j) = min(sim2.metric(1:j));
end

for k = 1:length(sim3.metric)
    best3(k) = min(sim3.metric(1:k));
end

stairs([1:length(sim1.metric)],best1,'b-')
hold on;
stairs([1:length(sim2.metric)],best2,'g-')
hold on;
stairs([1:length(sim3.metric)],best3,'r-')
legend('low', 'best', 'high');
