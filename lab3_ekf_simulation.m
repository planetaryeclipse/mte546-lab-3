%% lab3_ekf_simulation.m
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

% ekf update rate (must ensure is at least Nyquist)
ekf_freq = sampling_freq;
T_ekf = 1/ekf_freq;

% process update matrix
A = [1 T_ekf ; 0 1];

% process noise covariance
alpha = 1/T_ekf;
expected_x_std = 0.5; % mm
Q = [1 alpha ; alpha alpha^2 ]*expected_x_std^2;

% measurement noise (make sure matches dimension of z)
R = eye(2)*sensor_noise_var;
R1 = eye(2)*0.02;

% for convention
f = @(x) A*x;
h = @(x) [h_sensor(x(1)) ; h_sensor(x(1))];

F = @(x) A;
% 2 measurements both made with the same type of sensor
H = @(x) [(-sensor_b*sensor_c + sensor_a*sensor_d)/(sensor_c*x(1)+sensor_d)^2 0 ; ...
           (-sensor_b*sensor_c + sensor_a*sensor_d)/(sensor_c*x(1)+sensor_d)^2 0 ];

%% simulated trajectories

tf = 20;
simulated_time = 0:T_ekf:tf;

w_rand = mvnrnd(zeros(2,1),Q,length(simulated_time))';
eta_rand = normrnd(0,sqrt(R1(1,1)),2,length(simulated_time)); % same var

% linear point-to-point
xf = [200 0]';
x0 = [170 0]';
x_lin_pp = (xf-x0)/tf .* simulated_time + x0 + w_rand;
z_lin_pp = h_sensor(x_lin_pp(1,:)) + eta_rand;

% nonlinear point-to-point
A_nl_pp = 0.6*(xf-x0);
f_nl_pp = 9;
x_nl_pp = x_lin_pp + A_nl_pp*sin(2*pi*f_nl_pp/tf*simulated_time);
z_nl_pp = h_sensor(x_nl_pp(1,:)) + eta_rand;

% random motion (note must keep far from zero due to asymptote of h)
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
xlabel("Time (s)")
ylabel("Distance (mm)")
title("Simulated Trajectories")
legend("Linear PP", "Nonlinear PP", "Random Motion")

%% simulation

t = simulated_time;

% linear point-to-point

x_init_lin_pp = [170 0]';
P_init_lin_pp = Q;

[x_ests_lin_pp,P_ests_lin_pp] = evaluate_ekf(x_init_lin_pp,P_init_lin_pp,z_lin_pp,f,h,F,H,Q,R);
[lin_pp_stdPos, lin_pp_stdVel] = Uncertainty(P_ests_lin_pp);


% nonlinear point-to-point

x_init_nl_pp = [170 0]';
P_init_nl_pp = Q;

[x_ests_nl_pp,P_ests_nl_pp] = evaluate_ekf(x_init_nl_pp,P_init_nl_pp,z_nl_pp,f,h,F,H,Q,R);
[nonlin_pp_stdPos, nonlin_pp_stdVel] = Uncertainty(P_ests_nl_pp);


% random motion

x_init_rand_motion = [300 0]';
P_init_rand_motion = Q;

[x_ests_rand_motion,P_ests_rand_motion] = evaluate_ekf(x_init_rand_motion,P_init_rand_motion,z_rand_motion,f,h,F,H,Q,R);
[rand_pp_stdPos, rand_pp_stdVel] = Uncertainty(P_ests_rand_motion);

figure

subplot(3, 2, 1);
plot(t,x_lin_pp(1,:));
hold on
grid on
plot(t,x_ests_lin_pp(1,:),'--');
xlabel("Time (s)")
ylabel("Distance (mm)")
title("Linear Point to Point EKF Simulation")
legend("Actual","Estimates")

subplot(3, 2, 2);
plot(t,lin_pp_stdPos);
hold on
grid on
plot(t,lin_pp_stdVel);
xlabel("Time (s)")
ylabel("Standard Deviation (mm^2)")
title("Linear Point to Point Uncertainty")
legend("Distance std","Velocity std")

subplot(3, 2, 3);
plot(t,x_nl_pp(1,:));
hold on
grid on
plot(t,x_ests_nl_pp(1,:),'--');
xlabel("Time (s)")
ylabel("Distance (mm)")
title("Non-Linear Point to Point EKF Simulation")
legend("Actual","Estimates")

subplot(3, 2, 4);
plot(t,nonlin_pp_stdPos);
hold on
grid on
plot(t,nonlin_pp_stdVel);
xlabel("Time (s)")
ylabel("Standard Deviation (mm^2)")
title("Non-Linear Point to Point Uncertainty")
legend("Distance std","Velocity std")

subplot(3, 2, 5);
plot(t,x_rand_motion(1,:));
hold on
grid on
plot(t,x_ests_rand_motion(1,:),'--');
xlabel("Time (s)")
ylabel("Distance (mm)")
title("Random Motion EKF Simulation")
legend("Actual","Estimates")

subplot(3, 2, 6);
plot(t,rand_pp_stdPos);
hold on
grid on
plot(t,rand_pp_stdVel);
xlabel("Time (s)")
ylabel("Standard Deviation (mm^2)")
title("Non-Linear Point to Point Uncertainty")
legend("Distance std","Velocity std")

sgtitle("EKF Noise Covariance Sensors 0.02");

function [stdPos, stdVel] = Uncertainty(P)
    [~, ~, m] = size(P);
    
    stdPos = zeros(1, m); % Uncertainty Position
    stdVel = zeros(1, m); % Uncertainty Velocity
    
    % Compute standard deviations
    for k = 1:m
        stdPos(k) = sqrt(P(1, 1, k));
        stdVel(k) = sqrt(P(2, 2, k));
    end
end

