clc
clear
close all
filename = 'C:\Users\MNandyala\OneDrive - Inside MD Anderson\Documents\GitHub\fmincon_Entropy_arbitraryFunction\runs81.xlsx';
data = readmatrix(filename);
%%% Variable Setup
t_end = 1800; % [s] End time of the simulation
deltat = 60; % [s] Time step in which H is constant
% dt = 1; % time-step for forward run simulation of pennes solver
Ntime = t_end/deltat; % [-] Number of time step

%%% Time varying Magnetic Field Amplitude Vector Creation
Hmin = 7957; % [A/m] Minimum value of Magnetic Field the device can generate
Hmax = 39788; % [A/m] Maximum value of Magnetic Field the device can generate
H_range = (Hmin:Hmin:Hmax)'; % Possible steps of Magnetic Field on the device

InitialGuess(:,1) = 0*ones(Ntime,1);
InitialGuess(:,2) = H_range(1)*ones(Ntime,1);
InitialGuess(:,3) = H_range(2)*ones(Ntime,1);
InitialGuess(:,4) = H_range(3)*ones(Ntime,1);
InitialGuess(:,5) = H_range(4)*ones(Ntime,1);
InitialGuess(:,6) = H_range(5)*ones(Ntime,1);

% indices = repmat([1 2 3 4 5],1,Ntime/5);
indices = [3 5 1 2 4 2 3 5 1 3 2 4 5 2 1 4 5 3 1 5 4 2 3 1 5 2 3 1 5 2];
InitialGuess(:,7) = H_range(indices); % 1. Random variation of H bounded by the minimum and maximum values

[H_var_steps,] = H_variation_between_H_AND_0(Ntime);
InitialGuess(:,8) = 7600*H_var_steps;   % Rob's group H variation Peak amplitude was 7600 [A/m] with 50% duty cycle (60s On and 60s Off)

for i = 1:size(InitialGuess,2)

    Htime = InitialGuess(:,i); % Constant H in time

    for j = 1%:size(data,1)
        Qm_t = 0;                   % [W/m3] Tumor Metabolic Heat
        k_t=data(j,1); % [W/m/K] Thermal Conductivity
        w_t = data(j,2); % [1/s] Blood perfusion rate
        rhocp_t = data(j,3); % [J/m^3/K] Density and Specific Heat Production
        Km = data(j,4); % [J/m^3] Anisotropy Constant
        Gain(j,i) = pennesmht_fe(k_t,w_t,rhocp_t,Km,Htime,deltat);
        % Gains_InitialGuess(j,1:size(data,2)) = data(j,:);
        % Gains_InitialGuess(j,i+size(data,2)) = tempqoi;
    end
end
%% deltaT plots
% for run1
All_runs_initial = Temperatures(1,1,:);
All_runs_final = Temperatures(1,end,:);
DT = All_runs_final-All_runs_initial;
close all
plot(DT(:))
xlabel('Uncertain Variabl Combinations')
ylabel(['Change in Temperature (', char(176), 'C)'], 'FontSize', 14)
title('Concentration = 4 mg/cm^3')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15)  % Set for the current axes
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman') % Set for all elements in the figure
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 15) % Set font size for all elements in the figure
%%
SAR = [5.5,11,16.5,22,27.5];
plot(H_range,SAR)
SAR2 = SAR *1000
Q = SAR2*0.0192

%%
% close all
% % Define different line styles
% lineStyles = {'-', '--', ':', '-.', '-', '--', ':', '-.'};
% 
% % Plot the data
% subplot(3,3,1)
% hold on
% for i = 1:size(InitialGuess, 2)
%     stairs(InitialGuess(:, i), lineStyles{i}, 'LineWidth', 2);  % Apply different line styles and LineWidth=2
% end
% hold off;
% xlabel('Time Steps')
% ylabel('Magnetic Field, [A/m]')
% title('Different Initial Guesses')
% legend('IG1, H = 0','IG2, H = 7957','IG3, H = 15917','IG4, H = 23871','IG5, H = 31828','IG6, H = 39785','IG7, Random H','IG8, Cyclic H',location = 'westoutside')
% legendArray = {'IG1, H = 0','IG2, H = 7957','IG3, H = 15917','IG4, H = 23871','IG5, H = 31828','IG6, H = 39785','IG7, Random H','IG8, Cyclic H'};
% for i = 1:size(InitialGuess,2)
% subplot(3,3,i+1)
% plot(Gains_InitialGuess(:,i),lineStyles{i},'LineWidth', 2)
% xlabel('Runs')
% ylabel('Gain')
% title(['Gain variation, for ',legendArray{i}])
% end
% set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman', 'FontSize', 15);
%% Rob's Group H variation
function [rob_group_H,steps] = H_variation_between_H_AND_0(Ntime)
% Generate x values with a frequency of 1, starting from 1
steps = 1:1:Ntime;
% Calculate y values using the modified rectangle_wave function
rob_group_H = rectangle_wave(steps);
% Function definitions
    function y = rectangle_wave(x)
        % Define the rectangle wave function with cycle time of 2
        % Starting from upper bound at x = 1
        y = double(mod(x-1, 2) < 1);
    end
end
