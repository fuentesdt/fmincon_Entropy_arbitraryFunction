function [rob_group_H,steps] = H_variation_between_H_AND_0(Ntime)
% Generate x values with a frequency of 1, starting from 1
steps = 1:1:Ntime;

% Calculate y values using the modified rectangle_wave function
rob_group_H = rectangle_wave(steps);
%% Plot H variation
% figure;
% stairs(Tsteps, rob_group_H, 'LineWidth', 2);  % Using stairs for a more accurate representation
% title('Rectangle Wave (Cycle Time = 2, Frequency = 1)');
% xlabel('Time Step Number');
% ylabel('H Variation');
% ylim([-0.1, 1.1]);  % Set y-axis limits for better visibility
% xlim([1, Ntime]);   % Set x-axis limits, starting from 1
% grid on;
% 
% % Improve the appearance
% set(gca, 'FontSize', 12);
% set(gcf, 'Position', [100, 100, 800, 500]);

% Function definitions
function y = rectangle_wave(x)
    % Define the rectangle wave function with cycle time of 2
    % Starting from upper bound at x = 1
    y = double(mod(x-1, 2) < 1);
end

end