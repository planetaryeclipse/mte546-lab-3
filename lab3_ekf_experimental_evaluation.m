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

sensor_coeffs = [0.0295, 43.1601, 0.0943, 0.6177];
sensor_noise_var = noise_var;

clear coeffs;
clear noise_var;

% sensor model (x is mm)
sensor_a = sensor_coeffs(1);
sensor_b = sensor_coeffs(2);
sensor_c = sensor_coeffs(3);
sensor_d = sensor_coeffs(4);
h_sensor = @(x) (sensor_a*x+sensor_b)./(sensor_c*x+sensor_d);

distMeasured = @(x) (sensor_b - x*sensor_d)./(sensor_c*x - sensor_a);

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
std_filt_param = 0.1;
std_filt_movmean_k = 10;

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
% plot the unfiltered vs filtered data
figure;
plot(stat_time1, stat_data1(1, :), 'b', 'DisplayName', 'Row 1 - Original'); % Blue solid line
hold on;
plot(stat_time1, stat_data1(2, :), 'r', 'DisplayName', 'Row 2 - Original'); % Red solid line
plot(stat_time1, filt_stat_data1(1, :), 'b--', 'LineWidth', 2, 'DisplayName', 'Row 1 - Filtered'); % Blue dashed line
plot(stat_time1, filt_stat_data1(2, :), 'r--', 'LineWidth', 2, 'DisplayName', 'Row 2 - Filtered'); % Red dashed line
grid on;
xlabel('Time'); % Label for x-axis
ylabel('Data Value'); % Label for y-axis
title('Original vs Filtered Data'); % Title of the plot

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

%figure
%plot(stat_time1,stat_data1)

%% run the ekf on all experimental data
P0 = Q;

[stat_x_ests1,stat_P_ests1] = evaluate_ekf(stat_data1_x0,P0,filt_stat_data1,f,h,F,H,Q,R);
%[stat_x_ests2,stat_P_ests2] = evaluate_ekf(stat_data2_x0,P0,stat_data2,f,h,F,H,Q,R);
[stat_stdPos, stat_stdVel] = Uncertainty(stat_P_ests1);

stat_x_measured = distMeasured(stat_data1);

cropped_stat_time = stat_time1(15:150);
cropped_stat_ests = stat_x_ests1(1,15:150);

% Stationary EKF Performance and Uncertainty
figure
subplot(2, 1, 1);
plot(cropped_stat_time, cropped_stat_ests) % 35cm
hold on
grid on
plot(stat_time1, stat_x_measured) % 35cm
% Add markers at the start and end
plot(cropped_stat_time(1), 350, 'r.', 'MarkerSize', 20, 'DisplayName', 'Start'); % Start marker
plot(cropped_stat_time(end), 350, 'b.', 'MarkerSize', 20, 'DisplayName', 'End'); % End marker
xlabel("Time (s)")
ylabel("Distance (mm)")
legend("Estimated Position", "Measured Position", "Measured Start", "Measured End", 'Location', 'northwest'); % Moved to top left
title("Stationary EKF Performance")
subplot(2, 1, 2);
plot(stat_time1, stat_stdPos);
hold on
grid on
plot(stat_time1, stat_stdVel);
xlabel("Time (s)")
ylabel("Standard Deviation (mm^2)")
title("Stationary Uncertainty")
legend("Distance std", "Velocity std", 'Location', 'northwest'); % Moved to top left

[smooth_p2p_x_ests1,smooth_p2p_P_ests1] = evaluate_ekf(smooth_p2p_data1_x0,P0,filt_smooth_p2p_data1,f,h,F,H,Q,R);
%[smooth_p2p_x_ests2,smooth_p2p_P_ests2] = evaluate_ekf(smooth_p2p_data2_x0,P0,smooth_p2p_data2,f,h,F,H,Q,R);
[lin_pp_stdPos, lin_pp_stdVel] = Uncertainty(smooth_p2p_P_ests1);

cropped_smooth_time = smooth_p2p_time1(15:225);
cropped_smooth_ests = smooth_p2p_x_ests1(1,15:225);

% Linear Point to Point EKF Performance and Uncertainty
figure
subplot(2, 1, 1);
plot(cropped_smooth_time, cropped_smooth_ests) %35cm - 50cm
hold on
grid on
plot(cropped_smooth_time(1), 350, 'r.', 'MarkerSize', 20, 'DisplayName', 'Start'); % Start marker
plot(cropped_smooth_time(end), 500, 'b.', 'MarkerSize', 20, 'DisplayName', 'End'); % End marker
xlabel("Time (s)")
ylabel("Distance (mm)")
legend("Estimated Position", "Measured Start", "Measured End", 'Location', 'northwest'); % Moved to top left
title("Linear Point to Point EKF Performance")
subplot(2, 1, 2);
plot(smooth_p2p_time1, lin_pp_stdPos);
hold on
grid on
plot(smooth_p2p_time1, lin_pp_stdVel);
xlabel("Time (s)")
ylabel("Standard Deviation (mm^2)")
title("Linear Point to Point Uncertainty")
legend("Distance std", "Velocity std", 'Location', 'northwest'); % Moved to top left

[smooth_nonlinear_p2p_x_ests1,smooth_nonlinear_p2p_P_ests1] = evaluate_ekf(smooth_nonlinear_p2p_data1_x0,P0,filt_smooth_nonlinear_p2p_data1,f,h,F,H,Q,R);
%[smooth_nonlinear_p2p_x_ests2,smooth_nonlinear_p2p_P_ests2] = evaluate_ekf(smooth_nonlinear_p2p_data2_x0,P0,smooth_nonlinear_p2p_data2,f,h,F,H,Q,R);
[nonlin_pp_stdPos, nonlin_pp_stdVel] = Uncertainty(smooth_nonlinear_p2p_P_ests1);

cropped_nonlinear_time = smooth_nonlinear_p2p_time1(15:225);
cropped_nonlinear_ests = smooth_nonlinear_p2p_x_ests1(1,15:225);

% Nonlinear Point to Point EKF Performance and Uncertainty
figure
subplot(2, 1, 1);
hold on
grid on
plot(cropped_nonlinear_time(1), 350, 'r.', 'MarkerSize', 20, 'DisplayName', 'Start'); % Start marker
plot(cropped_nonlinear_time(end), 500, 'b.', 'MarkerSize', 20, 'DisplayName', 'End'); % End marker
plot(cropped_nonlinear_time, cropped_nonlinear_ests); %35cm - 50cm
xlabel("Time (s)")
ylabel("Distance (mm)")
legend("Estimated Position", "Measured Start", "Measured End", 'Location', 'northwest'); % Moved to top left
title("Nonlinear Point to Point EKF Performance")
subplot(2, 1, 2);
plot(smooth_nonlinear_p2p_time1, nonlin_pp_stdPos);
hold on
grid on
plot(smooth_nonlinear_p2p_time1, nonlin_pp_stdVel);
xlabel("Time (s)")
ylabel("Standard Deviation (mm^2)")
title("Non-Linear Point to Point Uncertainty")
legend("Distance std", "Velocity std", 'Location', 'northwest'); % Moved to top left


[random_x_ests1,random_P_ests1] = evaluate_ekf(random_data1_x0,P0,filt_random_data1,f,h,F,H,Q,R);
%[random_x_ests2,random_P_ests2] = evaluate_ekf(random_data2_x0,P0,filt_random_data2,f,h,F,H,Q,R);
[rand_stdPos, rand_stdVel] = Uncertainty(random_P_ests1);

cropped_random_time = random_time1(30:450);
cropped_random_ests = random_x_ests1(1,30:450);

% Random Motion EKF Performance and Uncertainty
figure
subplot(2, 1, 1);
plot(cropped_random_time, cropped_random_ests); % 60cm - 100cm
hold on
grid on
plot(cropped_random_time(1), 600, 'r.', 'MarkerSize', 20, 'DisplayName', 'Start'); % Start marker
plot(cropped_random_time(end), 1000, 'b.', 'MarkerSize', 20, 'DisplayName', 'End'); % End marker
xlabel("Time (s)")
ylabel("Distance (mm)")
legend("Estimated Position", "Measured Start", "Measured End", 'Location', 'northwest'); % Moved to top left
title("Random Motion EKF Performance")
subplot(2, 1, 2);
plot(random_time1, rand_stdPos);
hold on
grid on
plot(random_time1, rand_stdVel);
xlabel("Time (s)")
ylabel("Standard Deviation (mm^2)")
title("Random Motion Uncertainty")
legend("Distance std", "Velocity std", 'Location', 'northwest'); % Moved to top left

[tilted_x_ests1,tilted_P_ests1] = evaluate_ekf(tilted_data1_x0,P0,filt_tilted_data1,f,h,F,H,Q,R);
%[tilted_x_ests2,tilted_P_ests2] = evaluate_ekf(tilted_data2_x0,P0,filt_tilted_data2,f,h,F,H,Q,R);
[tilt_stdPos, tilt_stdVel] = Uncertainty(tilted_P_ests1);

cropped_tilted_time = tilted_time1(15:120);
cropped_tilted_ests = tilted_x_ests1(1,15:120);

% Tilted Point to Point EKF Performance and Uncertainty
figure
subplot(2, 1, 1);
plot(cropped_tilted_time, cropped_tilted_ests); % 35cm - 52.5cm
hold on
grid on
plot(cropped_tilted_time(1), 350, 'r.', 'MarkerSize', 20, 'DisplayName', 'Start'); % Start marker
plot(cropped_tilted_time(end), 525, 'b.', 'MarkerSize', 20, 'DisplayName', 'End'); % End marker
xlabel("Time (s)")
ylabel("Distance (mm)")
legend("Estimated Position", "Measured Start", "Measured End", 'Location', 'northwest'); % Moved to top left
title("Tilted Point to Point EKF Performance")
subplot(2, 1, 2);
plot(tilted_time1, tilt_stdPos);
hold on
grid on
plot(tilted_time1, tilt_stdVel);
xlabel("Time (s)")
ylabel("Standard Deviation (mm^2)")
title("Tilted Point to Point Uncertainty")
legend("Distance std", "Velocity std", 'Location', 'northwest'); % Moved to top left

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