%% lab3_ekf_formulation.m
% author: Samuel Street
% date: 2025-03-05

clear
clc
close all

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
sensor_a = sensor_coeffs(1);
sensor_b = sensor_coeffs(2);
sensor_c = sensor_coeffs(3);
sensor_d = sensor_coeffs(4);
h_sensor = @(x) (sensor_a*x+sensor_b)./(sensor_c*x+sensor_d);

%% motion model

% ekf update rate (must ensure is under Nyquist frequency)
nyquist_freq = sampling_freq/2;
T_ekf = 1/nyquist_freq;

% process update matrix
A = [1 T_ekf 1/2*T_ekf ; 0 1 T_ekf ; 0 0 1];

% process noise covariance
alpha = 1/T_ekf;
expected_x_std = 0.1; % mm
Q = [1 alpha 1/2*alpha ; alpha alpha^2 1/2*alpha^2 ; ...
    1/2*alpha 1/2*alpha^2 1/4*alpha^4]*expected_x_std^2;

Q = diag([0.2 0.4 1.0])*expected_x_std^2;

% measurement noise
R = sensor_noise_var;

% for convention
f = @(x) A*x;
h = @(x) h_sensor(x(1));

F = @(x) A;
H = @(x) [(-sensor_b*sensor_c + sensor_a*sensor_d)/(sensor_c*x(1)+sensor_d)^2 0 0];

%% simulated trajectories

tf = 20;
simulated_time = 0:T_ekf:tf;

w_rand = mvnrnd(zeros(3,1),Q,length(simulated_time))';
eta_rand = normrnd(0,sqrt(R),1,length(simulated_time));

% note that the velocity and acceleration are irrelevant
xf = [200 0 0]';
x0 = [150 0 0]';

% linear point-to-point
x_lin_pp = (xf-x0)/tf .* simulated_time + x0 + w_rand;
z_lin_pp = h_sensor(x_lin_pp(1,:)) + eta_rand;

% nonlinear point-to-point
A_nl_pp = 0.6*(xf-x0);
f_nl_pp = 9;
x_nl_pp = x_lin_pp + A_nl_pp*sin(2*pi*f_nl_pp/tf*simulated_time);
z_nl_pp = h_sensor(x_nl_pp(1,:)) + eta_rand;

% random motion
N = 3;
x_rand_motion = x0 + w_rand;
for k=1:N
    Ak = rand*(x0/2);
    Bk = rand*(x0/2);

    x_rand_motion = x_rand_motion ...
        + Ak*cos(2*pi*k/tf*simulated_time) ...
        + Bk*sin(2*pi*k/tf*simulated_time);
end
z_rand_motion = h_sensor(x_rand_motion(1,:)) + eta_rand;

%% plot the random trajectories
figure
plot(simulated_time,x_lin_pp(1,:));
hold on
grid on
plot(simulated_time,x_nl_pp(1,:));
plot(simulated_time,x_rand_motion(1,:));

legend("Linear PP", "Nonlinear PP", "Random Motion")

%% simulation

t = simulated_time;

% linear point-to-point

x_init_lin_pp = x0;
P_init_lin_pp = Q;

[x_ests_lin_pp,P_ests_lin_pp] = evaluate_ekf(x_init_lin_pp,P_init_lin_pp,z_lin_pp,f,h,F,H,Q,R);

figure
plot(t,x_lin_pp(1,:));
hold on
grid on
plot(t,x_ests_lin_pp(1,:),'--');

legend("Actual","Estimates")

% nonlinear point-to-point

x_init_nl_pp = x0;
P_init_nl_pp = Q;

[x_ests_nl_pp,P_ests_nl_pp] = evaluate_ekf(x_init_nl_pp,P_init_nl_pp,z_nl_pp,f,h,F,H,Q,R);

figure
plot(t,x_nl_pp(1,:));
hold on
grid on
plot(t,x_ests_nl_pp(1,:),'--');

legend("Actual","Estimates")

% random motion

x_init_rand_motion = x0;
P_init_rand_motion = Q;

[x_ests_rand_motion,P_ests_rand_motion] = evaluate_ekf(x_init_rand_motion,P_init_rand_motion,z_rand_motion,f,h,F,H,Q,R);

figure
plot(t,x_rand_motion(1,:));
hold on
grid on
plot(t,x_ests_rand_motion(1,:),'--');

legend("Actual","Estimates")

%% helper functions

function [x_ests,P_ests] = evaluate_ekf(x0,P0,z_actual,f,h,F,H,Q,R)

n = size(x0,1);
m = size(z_actual,2);

x_ests = zeros(n,m);
P_ests = zeros(n,n,m);

x_ests(:,1) = x0;
P_ests(:,:,1) = P0;

for k=2:m
    x_prev = x_ests(:,k-1);
    P_prev = P_ests(:,:,k-1);

    xk_predict = f(x_prev);
    zk_predict = h(x_prev);

    Fk = F(x_prev);
    Pk_predict = Fk*P_prev*Fk' + Q;

    Hk = H(xk_predict);
    Kk = Pk_predict*Hk'*(Hk*Pk_predict*Hk'+R)^-1;

    x_ests(:,k) = xk_predict+Kk*(z_actual(:,k)-zk_predict);
    P_ests(:,:,k) = (eye(n)-Kk*Hk)*Pk_predict;
end
end
