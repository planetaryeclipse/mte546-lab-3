std_filt_param = 0.1;
std_filt_movmean_k = 10;

filt_stat_data1 = std_clip(stat_data1,std_filt_param,std_filt_movmean_k);

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

function [y] = std_clip(data,std_dev,k)
mu = movmean(data,k,1);
std = sqrt(movvar(data,k,1));
y = max(min(mu+std*std_dev,data),mu-std*std_dev);
end