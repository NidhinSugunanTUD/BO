 function [cost] = runMPC(weights, param, dim, model, realmodel, do_plot)
if nargin < 6
    do_plot = 0;
end

% Extract weights
Qalpha = weights(1);
Qbeta = weights(2);
Ralpha = weights(3);
Rbeta = weights(4);

% Import CasADi
addpath('C:\Users\Nidhin\Downloads\Software\casadi-3.6.3-windows64-matlab2018b')
import casadi.*

%% Initialize Casadi
x = SX.sym('x', dim.nx); % state
u = SX.sym('u', dim.nu); % input
X = SX.sym('X', dim.nx, dim.N+1); % state trajectory
U = SX.sym('U', dim.nu, dim.N); % input trajectory
P = SX.sym('P', dim.nx + dim.nx); % parameters
X0 = P(1:6); % initial state
Ref = P(7:12); % state reference

%% Build Optimization Problem fields
% Dynamics
realdynamics = Function('dynamics', {x, u}, {realmodel.Ad*x + realmodel.Bd*u});

dynamics = Function('dynamics', {x, u}, {model.Ad*x + model.Bd*u});

% Symbolic solution for Single Shooting
X(:,1) = X0;
for k = 1:dim.N
    state = X(:,k);
    input = U(:,k);
    state_next  = dynamics(state,input);
    X(:,k+1) = state_next;
end
ssSol=Function('ssSol', {U, P}, {X});

% Objective function
Q = diag([Qalpha, 0, 0, Qbeta, 0, 0]);  % State weights
R = diag([Ralpha, Rbeta]);  % Input weights

objective = 0;
for k = 1:dim.N
    state = X(:,k);
    input = U(:,k);
    objective = objective + ...
        (state - Ref)' * Q * (state - Ref) + ...
        input' * R * input;
end

% Constraints
constraints = [];
for k = 1:dim.N+1
    constraints = [constraints ; X(1,k)];
    constraints = [constraints ; X(2,k)];
    constraints = [constraints ; X(3,k)];
    constraints = [constraints ; X(4,k)];
    constraints = [constraints ; X(5,k)];
    constraints = [constraints ; X(6,k)];
end

% Bounds
bnds = struct;

% state constraints
bnds.lbg(1:6:6*dim.N+6) = -pi;
bnds.lbg(2:6:6*dim.N+6) = -2*pi;
bnds.lbg(3:6:6*dim.N+6) = -50;
bnds.lbg(4:6:6*dim.N+6) = -pi;
bnds.lbg(5:6:6*dim.N+6) = -2*pi;
bnds.lbg(6:6:6*dim.N+6) = -50;

bnds.ubg(1:6:6*dim.N+6) = pi;
bnds.ubg(2:6:6*dim.N+6) = 2*pi;
bnds.ubg(3:6:6*dim.N+6) = 50;
bnds.ubg(4:6:6*dim.N+6) = pi;
bnds.ubg(5:6:6*dim.N+6) = 2*pi;
bnds.ubg(6:6:6*dim.N+6) = 50;

% input constraints
bnds.lbx(1:2:2*dim.N-1,1) = -10;
bnds.lbx(2:2:2*dim.N,1)   = -10;
bnds.ubx(1:2:2*dim.N-1,1) = 10;
bnds.ubx(2:2:2*dim.N,1)   = 10;

%%  Optimization problem
% initial guess
uGuess = zeros(dim.N,2);  % two control inputs

% fields
OPT_variables = reshape(U, 2*dim.N, 1); % single shooting
nlp_prob = struct;
nlp_prob.f = objective;
nlp_prob.x = OPT_variables;
nlp_prob.g = constraints;
nlp_prob.p = P;

% options
opts = struct;
opts.ipopt.max_iter = 200;
opts.ipopt.print_level = 0;
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-10;
opts.ipopt.acceptable_obj_change_tol = 1e-8;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

%% Initialize simulation
sim.xhat = zeros(dim.nx, dim.Nsim);
sim.x = zeros(dim.nx, dim.Nsim);
sim.u = zeros(dim.nu, dim.Nsim);

% with horizon
sim.xN = zeros(dim.nx, dim.N, dim.Nsim);
sim.uN = zeros(dim.nu, dim.N, dim.Nsim);

%% Simulation
x0 = [0; 0; 0; 0; 0; 0];     % initial condition
xref = [0.2; 0; 0; 0.3; 0; 0]; % Reference

sim.x(:, 1) = x0;
sim.xhat(:, 1) = x0;

for k = 1:dim.Nsim
    % set x0 and the arguments
    xk = sim.x(:, k);
    bnds.p   = [xk; xref]; % parameters
    bnds.x0 = reshape(uGuess', 2*dim.N, 1);

    % solve the MPC problem at step k
    sol = solver('x0', bnds.x0, 'lbx', bnds.lbx, 'ubx', bnds.ubx, ...
        'lbg', bnds.lbg, 'ubg', bnds.ubg,'p',bnds.p);

    uN = reshape(full(sol.x)', 2, dim.N)';
    xN = ssSol(uN',bnds.p); % compute optimal solution
    xN = full(xN);

    % extract first input
    u0 = uN(1,:)';

    % save the results
    sim.u(:, k) = u0;
    sim.xN(:,:,k) = xN(:, 1:dim.N);
    sim.uN(:,:,k) = uN';

    % update state
    if k < dim.Nsim
        % real model input and state
        xnext = realdynamics(xk, u0);
        sim.x(:, k+1) = full(xnext);

        % model input and state
        xhatnext = dynamics(xk, u0);
        sim.xhat(:, k+1) = full(xhatnext);
    end

end


%% Compute optimization metrics

% Objective fuction from paper...
    % performce oriented model learning for data driven MPC
        % equation 19
% cost = 10*ov^2 + 2* (sim.x(1, end) - xref(1))^2 + 2 * (sim.x(4, end) - xref(4))^2;

cost = log(mean(0.1*abs(xref(1)-sim.x(1,:)) + 0.1*abs(xref(4)-sim.x(4,:))) + 1) + ...
        log(0.1*(abs(xref(2)-sim.x(2, end))+ 0.1*(abs(xref(5)-sim.x(5, end)))) + 1);

%% Plot
if do_plot ~= 0
    t = [0:dim.Nsim] * param.Ts;
    alphaPlot = reshape(sim.x(1,:),1,[]);
    betaPlot = reshape(sim.x(4,:),1,[]);
    alphaDotPlot = reshape(sim.x(2,:),1,[]);
    betaDotPlot = reshape(sim.x(5,:),1,[]);

    f = figure(101);
    f.Position = [300 150 900 350];
    subplot(3,1,1)
    plot(t(1:end-1), alphaPlot, 'b', 'LineWidth', 1.5); hold on;
    plot(t(1:end-1), betaPlot, 'r', 'LineWidth', 1.5);
    plot(t(1:end-1), xref(1)*ones(1, dim.Nsim), 'b--', 'LineWidth', 2);
    plot(t(1:end-1), xref(4)*ones(1, dim.Nsim), 'r--', 'LineWidth', 2);
    
    legend('\alpha', '\beta', '\alpha_{ref}', '\beta_{ref}')
    xlabel('Time (seconds)')
    ylabel('angle (rad)')
    title('State Trajectories')
    box on


    subplot(3,1,2)
    plot(t(1:end-1), alphaDotPlot, 'b', 'LineWidth', 1.5); hold on;
    plot(t(1:end-1), betaDotPlot, 'r', 'LineWidth', 1.5);
    plot(t(1:end-1), zeros(1, dim.Nsim), 'k--', 'LineWidth', 1);
    
    legend({'$\dot{\alpha}$', '$\dot{\beta}$', 'Stationary'}, 'Interpreter', 'latex')
    xlabel('Time (seconds)')
    ylabel('angular velocity (rad/s)')
    title('State Trajectories')
    box on

    subplot(3,1,3)
    stairs(t(1:end-1), sim.u(1,:),'b','linewidth',1.5); hold on;
    stairs(t(1:end-1), sim.u(2,:),'r','linewidth',1.5);
    
    legend('\Deltap_\alpha', '\Deltap_\beta')
    ylabel('pressure (bar)')
    xlabel('Time (seconds)')
    title('Input Trajectories')
    box on
end
end