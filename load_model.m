% - the model is taken from
% Offset-free Model Predictive Control:
% A Ball Catching Application with a Spherical Soft Robotic Arm
% Yaohui Huang, Matthias Hofer and Raffaello D’Andrea

% general parameters
param.Ts = 0.02; %[s]

% model parameters
param.m = 4; %[kg]
param.radius = 10; %[cm]
param.g = 9.81; %[m/s^2]
param.k_alpha = 1.8; %[N*m/rad]
param.k_beta = 0.6;  %[N*m/rad]
param.d_alpha = 0.07; %[N*m*s/rad]
param.d_beta = 0.04;  %[N*m*s/rad]
param.h_alpha = 5.6; %[rad/m*bar]
param.h_beta = 3.7;  %[rad/m*bar]
param.tau_alpha = 0.003; %[s]
param.tau_beta = 0.004;  %[s]
param.c_alpha = 1e-5; %[m^2/(Nm)^2]
param.c_beta = 1e-6;  %[m^2/(Nm)^2]

%% Model:
% Continuous 
model.Ac = [     0,              1,                 0,                   0,              0,              0;
           -param.k_alpha, -param.d_alpha, -param.h_alpha,              0,              0,              0;
                0,          param.c_alpha, -1/param.tau_alpha,          0,              0,              0;
                0,              0,                 0,                   0,              1,              0;
                0,              0,                 0,           -param.k_beta,    -param.d_beta,   -param.h_beta;
                0,              0,                 0,                   0,         param.c_alpha, -1/param.tau_alpha];
            
model.Bc = [     0,            0;
                0,            0;
           1/param.tau_alpha, 0;
                0,            0;
                0,            0;
                0,  1/param.tau_beta];
            
model.Cc = [1, 0, 0, 0, 0, 0
            0, 0, 0, 1, 0, 0];

model.sysc = ss(model.Ac, model.Bc, model.Cc, []);

% Discrete
model.sysd = c2d(model.sysc, param.Ts, 'zoh');

model.Ad = model.sysd.A;
model.Bd = model.sysd.B;
model.Cd = model.sysd.C;

%% Real Model
realparam = param;
realparam.m = realparam.m * 1.23; %[kg]
realparam.radius = realparam.radius * 1.19; %[cm]
realparam.k_alpha = realparam.k_alpha * 0.83; %[N*m/rad]
realparam.k_beta = realparam.k_beta * 0.88;  %[N*m/rad]
realparam.d_alpha = realparam.d_alpha * 1.15; %[N*m*s/rad]
realparam.d_beta = realparam.d_beta * 0.88;  %[N*m*s/rad]
realparam.h_alpha = realparam.h_alpha * 1.17; %[rad/m*bar]
realparam.h_beta = realparam.h_beta * 1.16;  %[rad/m*bar]
realparam.tau_alpha = realparam.tau_alpha * 1.25; %[s]
realparam.tau_beta = realparam.tau_beta * 0.88;  %[s]
realparam.c_alpha = realparam.c_alpha * 1.16; %[m^2/(Nm)^2]
realparam.c_beta = realparam.c_beta * 1.17;  %[m^2/(Nm)^2]

% Continuous 
realmodel.Ac = [     0,              1,                 0,                      0,              0,              0;
       -realparam.k_alpha, -realparam.d_alpha, -realparam.h_alpha,              0,              0,              0;
                0.1,          realparam.c_alpha, -1/realparam.tau_alpha,          0,              0,              0;
                0,              0,                      0,                      0,              1,              0;
                0,              0,                      0,           -realparam.k_beta,    -realparam.d_beta,   -realparam.h_beta;
                0,              0,                      0,                      0,         realparam.c_alpha, -1/realparam.tau_alpha];
    

% add noise
randuncert = [0.0250   -0.0446    0.0264    0.0425   -0.0535   -0.0134;
   -0.1108   -0.0209    0.0705    0.0774    0.0686    0.0568;
   -0.0631   -0.0005   -0.0120   -0.0368   -0.0103   -0.0382;
   -0.0937   -0.0874   -0.0091    0.0336   -0.0315    0.0436;
   -0.0099   -0.0422    0.0075   -0.0885   -0.0210    0.0120;
   -0.0338    0.0347   -0.0641    0.0515   -0.0102   -0.0247];
realmodel.Ac = realmodel.Ac + randuncert;

realmodel.Bc = [    0,                  0;
                    0,                  0;
               1/realparam.tau_alpha,   0;
                    0,                  0;
                    0,                  0;
                    0,        1/realparam.tau_beta];


            
realmodel.Cc = [1, 0, 0, 0, 0, 0
                0, 0, 0, 1, 0, 0];

realmodel.sysc = ss(realmodel.Ac, realmodel.Bc, realmodel.Cc, []);

% Discrete
realmodel.sysd = c2d(realmodel.sysc, param.Ts, 'zoh');

realmodel.Ad = realmodel.sysd.A;
realmodel.Bd = realmodel.sysd.B;
realmodel.Cd = realmodel.sysd.C;

%% Analysis

% % Pole-zero mapping
% figure(1)
% subplot(1,2,1)
% pzmap(model.sysd);
% title('Pole-zero map of model')
% subplot(1,2,2)
% pzmap(realmodel.sysd);
% title('Pole-zero map of process')
% 
% % Step response
% stepmodel = step(model.sysc);
% realstepmodel = step(realmodel.sysc);
% stepmode_11 = stepmodel(:,1,1);
% stepmode_22 = stepmodel(:,2,2);
% realstepmode_11 = realstepmodel(:,1,1);
% realstepmode_22 = realstepmodel(:,2,2);
% 
% figure(2)
% subplot(2,2,1)
% plot(stepmode_11)
% title('Model''s step response[1:1]')
% subplot(2,2,2)
% plot(stepmode_22)
% title('Model''s step response[2:2]')
% subplot(2,2,3)
% plot(realstepmode_11)
% title('Process''s step response[1:1]')
% subplot(2,2,4)
% plot(realstepmode_22)
% title('Process''s step response[2:2]')
% sys_info_11 = stepinfo(stepmode_11)
% sys_info_22 = stepinfo(stepmode_22)
% 
% % Check if poles are inside the unit circle
% poles = pole(model.sysd);
% realpoles = pole(realmodel.sysd);
% arePolesInside = all(abs(poles) < 1);
% if arePolesInside
%     disp('All poles are inside the unit circle. The model is stable.');
% else
%     disp('Some poles are outside the unit circle. The model is not stable.');
% end
% 
% areRealPolesInside = all(abs(realpoles) < 1);
% if areRealPolesInside
%     disp('All poles are inside the unit circle. The process is stable.');
% else
%     disp('Some poles are outside the unit circle. The process is not stable.');
% end
% 
% % Check controllability
% Co_model = ctrb(model.Ad,model.Bd);
% uncontrollable_states_model = length(model.Ad) - rank(Co_model) % Determine the number of uncontrollable states
% Co_process = ctrb(realmodel.Ad,realmodel.Bd);
% uncontrollable_states_process = length(realmodel.Ad) - rank(Co_process)


% % RGA
% modeltf = tf(model.sysc);
% realmodeltf = tf(realmodel.sysc);
% 
% omega = [0 0.4*2*pi];
% for i=1:length(omega)
%     RGA(:,:,i) = evalfr( modeltf, omega(i) ).* pinv( evalfr(modeltf,omega(i)) ).';
%     realRGA(:,:,i) = evalfr( realmodeltf, omega(i) ).* pinv( evalfr(realmodeltf,omega(i)) ).';
% end

%% Dimensions
dim.N = 20; % Prediction horizon
dim.nx = size(model.Ad, 1);
dim.nu = size(model.Bd, 2);
dim.Nsim = 30; % Simulation time