clc
clear
close all

% %%
% SpaceData = 'SpatialProfiles.xlsx';
% SpatialProfiles = readmatrix(SpaceData);
% TemporalData = 'TemporalProfiles.xlsx';
% TemporalProfiles = readmatrix(TemporalData);
%%
%%% Variable Setup
t_end = 1800; % [s] End time of the simulation
deltat = 60; % [s] Delta
Ntime = t_end/deltat; % [-] Number of step in H variation

%%% Mean uncertain variables
k_t=.492; % [W/m/K] Thermal Conductivity
w_t = 0.00682; % [1/s] Blood perfusion rate
rhocp_t = 3589229; % [J/m^3/K] Density and Specific Heat Production
Km = 31500; % [J/m^3] Anisotropy Constant %% A variable to be recoved from the experiments/data

%%% Time varying Magnetic Field Amplitude Vector Creation
Hmin = 8000; % [A/m] Minimum value of Magnetic Field the device can generate
Hmax = 40000; % [A/m] Maximum value of Magnetic Field the device can generate
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
    close all
    Htime = InitialGuess(:,i);
    tic
    Gain(i) = pennesmht_fe3o4(k_t,w_t,rhocp_t,Km,Htime,deltat,i);
    toc
end
% close all
% tic
% pennesmht_fe(k_t,w_t,rhocp_t,Km,Htime,deltat)
% toc

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