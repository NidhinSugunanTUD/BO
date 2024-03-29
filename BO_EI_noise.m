%% Bayesian Optimization 
close all; clc; clear all;

% Load 6 state LTI model
run load_model;
run load_noise.m;

%% Choose noise
% Choose noise between:
%   - zero_noise
%   - low_noise
%   - medium_noise
%   - high_noise
%
noise = high_noise;
Niter = 5; % should be dim.nnoise

Nstart = 1;     % Initial no. of observations
Nobs = 50;      % No. of more observations to perform

all_metrics = zeros(Niter, Nstart+Nobs);
for in = 1 : Niter
    noisei = noise(:,:,in);

    %% Setup BO
    % Random seed
    rng(1);
    warning('off','all');
    
    % Expected Improvement (EI)
    xi = 0.02;  % Exploration-exploitation parameter 0.4
    % MPC function
    objective_function = @(weights) runMPC(weights,param,dim,model,realmodel);
    
    % Search range for Optimization variables 
    Wrange = [200 5000; % Qalpha
              200 5000; % Qbeta
              0.005 1;    % Ralpha
              0.005 1];   % Rbeta
    
    [w1, w2, w3, w4] ...
            = ndgrid(linspace(Wrange(1,1),Wrange(1,2), 31),...
                      linspace(Wrange(2,1),Wrange(2,2), 31),...
                      linspace(Wrange(3,1),Wrange(3,2), 31),...
                      linspace(Wrange(4,1),Wrange(4,2), 31));
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
    % Intialize simulation arrays
    sim.W = zeros(4, Nstart+Nobs);
    sim.metric = zeros(1, Nstart+Nobs);
    
    % Starting point
    W0 = [300; 300; 1; 1]; %[800; 800; 1; 1]; 
    
    sim.W(:, 1) = W0;
    sim.metric(1) = objective_function(W0);
    
    %% Run BO
    hw = waitbar(0,'Running BO...');
    
    for j = 1:Nobs
        
        % Gaussian process model
        Wj = sim.W(:, 1:j);
        metricj = sim.metric(1:j);
        surrogate_function = fitrgp(Wj', metricj', ...
                                    'KernelFunction','squaredexponential');
        
        [mu, sigma] = predict(surrogate_function, W_grid); % mean and standard deviation
        
        % Expected Improvement
        imp = distance_to_opt(mu, metricj);
        Z = (sigma ~= 0).* imp ./ sigma;
        Phi = @(x) normcdf(x); % Normal distribution CDF
        phi = @(x) normpdf(x); % Normal distribution PDF
        EI = imp.* Phi(Z) +  sigma .* phi(Z); 
        
        % Take the point that gets the maximum EI 
        posEI = find(EI == max(EI));
    
        % If more points with the same EI, take a random one
        random_pos = randi(length(posEI));
        xEI = W_grid(posEI(random_pos), :);
        
        % Save next points
        % xEI_noise = xEI' ./ noise(:,Nobs,Niter);
        xEI_noise = xEI' ./ noisei(:,Nobs);
        % Save next points
        sim.W(:, j+1) = xEI_noise';
        sim.metric(j+1) = objective_function(xEI_noise);
            
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
    figure(1 + 1000*i)
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
    h = xlabel('Iteration index', 'FontSize', 14);  % label x axis
    set(h,'Interpreter', 'Latex');  % label
    % xticks([0, 25, 50]);
    % xlim([0, 52]);
    h = ylabel('Performance cost', 'FontSize', 14); % label y axis
    set(h,'Interpreter', 'Latex');  % label
    % yticks([0, 0.065, 0.15]);
    % ylim([0, 0.15]);
    % legend('Current point', 'Current best point', 'Overall best point');
    grid('on');

    all_metrics(in, :) = sim.metric;
   
end