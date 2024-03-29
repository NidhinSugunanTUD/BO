%% Bayesian Optimization 
close all; clc; clear all;

% Load 6 state LTI model
 run load_model;

%% Setup BO
% Random seed
rng(502460);
warning('off','all');

% Probability of Improvement (PI)
xi = 0;  % Exploration-exploitation parameter 

% MPC function
objective_function = @(weights) runMPC(weights,param,dim,model,realmodel);

% Search range for Optimization variables 
Wrange = [200 5000; % Qalpha
          200 5000; % Qbeta
          0.005 1;    % Ralpha
          0.005 1];   % Rbeta

[w1, w2, w3, w4] ...
        = ndgrid(linspace(Wrange(1,1),Wrange(1,2), 241),...
                  linspace(Wrange(2,1),Wrange(2,2), 241),...
                  linspace(Wrange(3,1),Wrange(3,2), 21),...
                  linspace(Wrange(4,1),Wrange(4,2), 21));
W_grid = [w1(:), w2(:), w3(:), w4(:)];

% Choose minimization (1) / maximization (0)
is_minimization = 1;
if is_minimization == 1
    % minimization
    distance_to_opt = @(mu, metric) min(metric) - mu - xi;
else 
    % maximization
    distance_to_opt = @(mu, metric) mu - max(metric) - xi;
end

%% Initialization
Nstart = 1;     % Initial no. of observationsx`
Nobs = 5;       % No. of more observations to perform

% Intialize simulation arrays
sim.W = zeros(4, Nstart+Nobs);
sim.metric = zeros(1, Nstart+Nobs);
sim.radcost = zeros(1, Nstart+Nobs);

% Starting point
W0 = [3760;     4580;       1;       0.95025];
sim.W(:, 1) = W0;
[cost0, radcost0] = objective_function(W0);
sim.metric(1) = cost0;
sim.radcost(1) = radcost0;

%% Run BO
hw = waitbar(0,'Running BO...');

for j = 1:Nobs
    
    % Gaussian process model
    Wj = sim.W(:, 1:j);
    metricj = sim.metric(1:j);
    surrogate_function = fitrgp(Wj', metricj', ...
                                'KernelFunction','squaredexponential');
    
    [mu, sigma] = predict(surrogate_function, W_grid);
    
    % Probability of Improvement
    imp = distance_to_opt(mu, metricj);
    Z = (sigma ~= 0).* imp ./ sigma;
    Phi = @(x) normcdf(x); % Normal distribution CDF
    PI = Phi(Z); 
    
    % Take the point that gets the maximum EI 
    posPI = find(PI == max(PI));

    % If more points with the same EI, take a random one
    random_pos = randi(length(posPI));
    xPI = W_grid(posPI(random_pos), :);
    
    % Save next points
    sim.W(:, j+1) = xPI';
    [cost, radcost] = objective_function(xPI);
    sim.metric(j+1) = cost;
    sim.radcost(j+1) = radcost;
        
    % Update progress bar
    waitbar(j/Nobs,hw);

end
close(hw)

%% Extract best point
if is_minimization
    str = 'Minimum';
    [ae,be] = min(mu);
    [ao,bo] = min(sim.metric); 
else
    str = 'Maximum';
    [ae,be] = max(mu); 
    [ao,bo] = max(sim.metric); 
end
fprintf('Bayesian Optimization\n');
fprintf('  %s (estimated): y(%.6f) = %.6f\n',str, W_grid(be),ae);
fprintf('  %s (observed) : y(%.6f) = %.6f\n',str, sim.W(bo),ao);

%% Display BO Table
names = {'Qalpha', 'Qbeta', 'Ralpha', 'Rbeta', 'J (cost)'};
botable = array2table([sim.W', sim.metric'], 'VariableNames', names)
botable_sorted = sortrows(botable, 'J (cost)')

%% Plot the winner
opt_weights = sim.W(:,bo);
runMPC(opt_weights, param, dim, model, realmodel, 1)

%% Final plot 
figure(1)
scatter([1:length(sim.metric)],sim.metric,'k*')
hold on;
best = zeros(1,length(sim.metric));

for i = 1:length(sim.metric)
    best(i) = min(sim.metric(1:i));
end

stairs([1:length(sim.metric)],best,'b-')
hold on
[minBEST, minIndex] = min(sim.metric);
scatter(minIndex, minBEST, 120, 'rs', 'filled')

% Labelling
h = xlabel('Iteration index'); %label x axis
set(h,'Interpreter', 'Latex'); %label
h = ylabel('Performance cost'); %label y axis
set(h,'Interpreter', 'Latex'); %label
legend('Current point', 'Current best point', 'Overall best point');
grid('on');