%% lab3_ekf_formulation.m
% author: Samuel Street
% date: 2025-03-05

clear
clc

%% loads the long sensor data

sampling_freq = 15;
T_sampling = 1/sampling_freq;

% both sensors used are long range IR sensors
load("long_coeffs.mat","coeffs")
load("long_noise.mat","noise_var")

sensor_coeffs = coeffs;
sensor_noise_var = noise_var;

clear coeffs;
clear noise_var;

% sensor model (x is mm)
h_sensor = @(x) (sensor_coeffs(1).*x + sensor_coeffs(2)) ./ (sensor_coeffs(3).*x + sensor_coeffs(4));

%% motion model

% ekf update rate (must ensure is under Nyquist frequency)
nyquist_freq = sampling_freq/2-1;
T_ekf = 1/nyquist_freq;

% process update matrix
A = [1 T_ekf 1/2*T_ekf^2 ; 0 1 T_ekf ; 0 0 1];

% process noise covariance
%alpha = 1e-2; % 1/T_ekf
alpha = 1/T_ekf;
expected_x_std = 2; % mm
Q = [1 alpha 1/2*alpha ; alpha alpha^2 1/2*alpha^2 ; ...
    1/2*alpha 1/2*alpha^2 1/4*alpha^4]*expected_x_std^2*1e-3;

% Q = diag([1e-3 1e-5 1e-6]);

% measurement matrix
H = [sensor_coeffs(1)/sensor_coeffs(3) 0 0];

% measurement noise
R = sensor_noise_var;

%% simulated trajectories

tf = 20;
simulated_time = 0:T_ekf:tf;

w_rand = mvnrnd(zeros(3,1),Q,length(simulated_time))';
eta_rand = normrnd(0,sqrt(R),1,length(simulated_time));

% linear point-to-point
xf = [200 0 0]';
x0 = [150 0 0]';
x_lin_pp = (xf-x0)/tf .* simulated_time + x0 + w_rand;
z_lin_pp = h_sensor(x_lin_pp(1,:)) + eta_rand;

% nonlinear point-to-point
A_nl_pp = 0.6*(xf-x0);
f_nl_pp = 9;
x_nl_pp = x_lin_pp + A_nl_pp*sin(2*pi*f_nl_pp/tf*simulated_time);
z_nl_pp = h_sensor(x_nl_pp(1,:)) + eta_rand;

% random motion
N = 10;
x_rand_motion = x0 + w_rand;
for k=1:N
    Ak = rand*(x0/2);
    Bk = rand*(x0/2);

    x_rand_motion = x_rand_motion ...
        + Ak*cos(2*pi*k/tf*simulated_time) ...
        + Bk*sin(2*pi*k/tf*simulated_time);
end
z_rand_motion = h_sensor(x_rand_motion(1,:)) + eta_rand;

% plot the random trajectories
figure
plot(simulated_time,x_lin_pp(1,:));
hold on
grid on
plot(simulated_time,x_nl_pp(1,:));
plot(simulated_time,x_rand_motion(1,:));

legend("Linear PP", "Nonlinear PP", "Random Motion")

%% set simulation parameters

% ekf model parameters

% -- already set --
% T_ekf
% sensor_coeffs
% Q
% R

x_init = zeros(3,1);
P_init = Q;

A_model = A;
F_linear = A;
H_linear = H;

%% simulation - linear point-to-point

t_actual = simulated_time;
z_actual_measure = z_lin_pp;
x_pos_actual = x_lin_pp(1,:);

simout_lin_pp = sim("lab3_ekf_simulation.slx","StopTime",num2str(tf));
x_pos_estimate_lin_pp = squeeze(simout_lin_pp.x_predict(1,:,:))';

figure
plot(simulated_time,x_lin_pp(1,:));
hold on
grid on
plot(simout_lin_pp.tout,x_pos_estimate_lin_pp,'--');
legend("Actual","Estimate")

%% simulation - nonlinear point-to-point

t_actual = simulated_time;
z_actual_measure = z_nl_pp;
x_pos_actual = x_nl_pp(1,:);

simout_nl_pp = sim("lab3_ekf_simulation.slx","StopTime",num2str(tf));
x_pos_estimate_nl_pp = squeeze(simout_nl_pp.x_predict(1,:,:))';

figure
plot(t_actual,x_nl_pp(1,:));
hold on
grid on
plot(t_actual,x_pos_estimate_nl_pp,'--');
legend("Actual","Estimate")

%% simulation - random motion

t_actual = simulated_time;
z_actual_measure = z_rand_motion;
x_pos_actual = x_rand_motion(1,:);

simout_rand_motion = sim("lab3_ekf_simulation.slx","StopTime",num2str(tf));
x_pos_estimate_rand_motion = squeeze(simout_rand_motion.x_predict(1,:,:))';

figure
plot(t_actual,x_rand_motion(1,:));
hold on
grid on
plot(t_actual,x_pos_estimate_rand_motion,'--');
legend("Actual","Estimate")