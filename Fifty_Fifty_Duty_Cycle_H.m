clc
clear
close all
Ntime = 30;
Tsteps = 1:1:Ntime;
[rob_group_H,~] = H_variation_between_H_AND_0(Ntime);
% Plot the rectangle wave
figure;
stairs(Tsteps, rob_group_H, 'LineWidth', 2);  % Using stairs for a more accurate representation
title('50% Duty Cycle');
xlabel('Time Step Number');
ylabel('H Variation');
ylim([-0.1, 1.1]);  % Set y-axis limits for better visibility
xlim([1, Ntime]);   % Set x-axis limits, starting from 1
grid on;

% Improve the appearance
set(gca, 'FontSize', 12);
set(gcf, 'Position', [100, 100, 800, 500]);


