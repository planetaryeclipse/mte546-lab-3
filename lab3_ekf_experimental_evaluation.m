%% lab3_ekf_experimental_evaluation.m
% author: Samuel Street
% date: 2025-03-06

clc
clear
close all

%% loads the long sensor data

sampling_freq = 15; % Hz
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

%% load all experimental data

[stat_data1,stat_time1,stat_data1_x0] = load_exp_data("collected_data/motion_stationary.mat");
[stat_data2,stat_time2,stat_data2_x0] = load_exp_data("collected_data/motion_stationary2.mat");

[smooth_p2p_data1,smooth_p2p_time1,smooth_p2p_data1_x0] = load_exp_data("collected_data/motion_smoothp2p.mat");
[smooth_p2p_data2,smooth_p2p_time2,smooth_p2p_data2_x0] = load_exp_data("collected_data/motion_smoothp2p2.mat");

[smooth_nonlinear_p2p_data1,smooth_nonlinear_p2p_time1,smooth_nonlinear_p2p_data1_x0] = load_exp_data("collected_data/motion_nonlinearp2p.mat");
[smooth_nonlinear_p2p_data2,smooth_nonlinear_p2p_time2,smooth_nonlinear_p2p_data2_x0] = load_exp_data("collected_data/motion_nonlinearp2p2.mat");

[random_data1,random_time1,random_data1_x0] = load_exp_data("collected_data/motion_random.mat");
[random_data2,random_time2,random_data2_x0] = load_exp_data("collected_data/motion_random2.mat");

[tilted_data1,tilted_time1,tilted_data1_x0] = load_exp_data("collected_data/motion_tilted.mat");
[tilted_data2,tilted_time2,tilted_data2_x0] = load_exp_data("collected_data/motion_tilted2.mat");

% filter the experimental data using std_filt
std_filt_param = 1;
std_filt_movmean_k = 30;

filt_stat_data1 = std_clip(stat_data1,std_filt_param,std_filt_movmean_k);
filt_stat_data2 = std_clip(stat_data2,std_filt_param,std_filt_movmean_k);

filt_smooth_p2p_data1 = std_clip(smooth_p2p_data1,std_filt_param,std_filt_movmean_k);
filt_smooth_p2p_data2 = std_clip(smooth_p2p_data2,std_filt_param,std_filt_movmean_k);

filt_smooth_nonlinear_p2p_data1 = std_clip(smooth_nonlinear_p2p_data1,std_filt_param,std_filt_movmean_k);
filt_smooth_nonlinear_p2p_data2 = std_clip(smooth_nonlinear_p2p_data2,std_filt_param,std_filt_movmean_k);

filt_random_data1 = std_clip(random_data1,std_filt_param,std_filt_movmean_k);
filt_random_data2 = std_clip(random_data2,std_filt_param,std_filt_movmean_k);

filt_tilted_data1 = std_clip(tilted_data1,std_filt_param,std_filt_movmean_k);
filt_tilted_data2 = std_clip(tilted_data2,std_filt_param,std_filt_movmean_k);

% plot the unfiltered vs filtered data
figure
plot(stat_time1,stat_data1(1,:))
hold on
grid on
plot(stat_time1,stat_data1(2,:))
plot(stat_time1,filt_stat_data1(1,:),"--","LineWidth",2)
plot(stat_time1,filt_stat_data1(2,:),"--","LineWidth",2)

legend("Unfiltered $z_1$","Unfiltered $z_2$","Filtered $z_1$","Filtered $z_2$","Interpreter","latex")

%% configure the ekf parameters
% note refer to lab3_ekf_simulation.m for more details on the motion model
% and how each parameter should be configured

ekf_freq = sampling_freq;
T_ekf = 1/ekf_freq;

A = [1 T_ekf ; 0 1];

alpha = 1/T_ekf;
expected_x_std = 0.5; % mm
Q = [1 alpha ; alpha alpha^2 ]*expected_x_std^2;
R = eye(2)*sensor_noise_var;

f = @(x) A*x;
h = @(x) [h_sensor(x(1)) ; h_sensor(x(1))];

F = @(x) A;
H = @(x) [(-sensor_b*sensor_c + sensor_a*sensor_d)/(sensor_c*x(1)+sensor_d)^2 0 ; ...
           (-sensor_b*sensor_c + sensor_a*sensor_d)/(sensor_c*x(1)+sensor_d)^2 0 ];

%% plot some of the measured observations

figure
plot(stat_time1,stat_data1)

%% run the ekf on all experimental data
P0 = Q;

[stat_x_ests1,stat_P_ests1] = evaluate_ekf(stat_data1_x0,P0,stat_data1,f,h,F,H,Q,R);
[stat_x_ests2,stat_P_ests2] = evaluate_ekf(stat_data2_x0,P0,stat_data2,f,h,F,H,Q,R);

figure
plot(stat_time1,stat_x_ests1(1,:))
hold on
grid on
plot(stat_time2,stat_x_ests2(1,:))
xlabel("Time (s)")
ylabel("Distance (m)")

[smooth_p2p_x_ests1,smooth_p2p_P_ests1] = evaluate_ekf(smooth_p2p_data1_x0,P0,smooth_p2p_data1,f,h,F,H,Q,R);
[smooth_p2p_x_ests2,smooth_p2p_P_ests2] = evaluate_ekf(smooth_p2p_data2_x0,P0,smooth_p2p_data2,f,h,F,H,Q,R);

figure
plot(smooth_p2p_time1,smooth_p2p_x_ests1(1,:))
hold on
grid on
plot(smooth_p2p_time2,smooth_p2p_x_ests2(1,:))
xlabel("Time (s)")
ylabel("Distance (m)")

[smooth_nonlinear_p2p_x_ests1,smooth_nonlinear_p2p_P_ests1] = evaluate_ekf(smooth_nonlinear_p2p_data1_x0,P0,smooth_nonlinear_p2p_data1,f,h,F,H,Q,R);
[smooth_nonlinear_p2p_x_ests2,smooth_nonlinear_p2p_P_ests2] = evaluate_ekf(smooth_nonlinear_p2p_data2_x0,P0,smooth_nonlinear_p2p_data2,f,h,F,H,Q,R);

figure
plot(smooth_nonlinear_p2p_time1,smooth_nonlinear_p2p_x_ests1(1,:));
hold on
grid on
plot(smooth_nonlinear_p2p_time2,smooth_nonlinear_p2p_x_ests2(1,:));
xlabel("Time (s)")
ylabel("Distance (m)")

[random_x_ests1,random_P_ests1] = evaluate_ekf(random_data1_x0,P0,random_data1,f,h,F,H,Q,R);
[random_x_ests2,random_P_ests2] = evaluate_ekf(random_data2_x0,P0,random_data2,f,h,F,H,Q,R);

figure
plot(random_time1,random_x_ests1(1,:));
hold on
grid on
plot(random_time2,random_x_ests2(1,:));
xlabel("Time (s)")
ylabel("Distance (m)")

[tilted_x_ests1,tilted_P_ests1] = evaluate_ekf(tilted_data1_x0,P0,tilted_data1,f,h,F,H,Q,R);
[tilted_x_ests2,tilted_P_ests2] = evaluate_ekf(tilted_data2_x0,P0,tilted_data2,f,h,F,H,Q,R);

figure
plot(tilted_time1,tilted_x_ests1(1,:));
hold on
grid on
plot(tilted_time2,tilted_x_ests2(1,:));
xlabel("Time (s)")
ylabel("Distance (m)s")

%% helpers

function [data,time,x0] = load_exp_data(filename)
% quick wrapper function for convenience and organization
load(filename,"data","time")

% ensure sequential data is provided columnwise
data = data';
time = time';

x0 = data(:,1); % can be used as the starting value of the ekf
end

function [y] = std_clip(data,std_dev,k)
mu = movmean(data,k,1);
std = sqrt(movvar(data,k,1));
y = max(min(mu+std*std_dev,data),mu-std*std_dev);
end