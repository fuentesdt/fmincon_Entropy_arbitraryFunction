clc
clear
close all
filename = 'C:\Users\MNandyala\OneDrive - Inside MD Anderson\Documents\GitHub\fmincon_Entropy_arbitraryFunction\runs81.xlsx';
data = readmatrix(filename);
%%
SpaceData = 'SpatialProfiles.xlsx';
SpatialProfiles = readmatrix(SpaceData);
TemporalData = 'TemporalProfiles.xlsx';
TemporalProfiles = readmatrix(TemporalData);
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
Hmin = 7957; % [A/m] Minimum value of Magnetic Field the device can generate
Hmax = 39788; % [A/m] Maximum value of Magnetic Field the device can generate
H_range = (Hmin:Hmin:Hmax)'; % Possible steps of Magnetic Field on the device
InitialGuess = H_range(5)*ones(Ntime,1);  IG =6;
Htime = InitialGuess; % Constant H in time

H = randi([Hmin, Hmax], 1, 500);


for i = 1: length(H)
    Htime = H(i)*ones(Ntime,1);
    gain(i) = pennesmht_fe(k_t,w_t,rhocp_t,Km,Htime,deltat);
end

close all
plot(H,gain,'o')
xlabel('K')
ylabel('Gain G(K,P) @ fixed P')


% tic
% pennesmht_fe3o4(k_t,w_t,rhocp_t,Km,Htime,deltat)
% toc
% close all
% tic
% pennesmht_fe(k_t,w_t,rhocp_t,Km,Htime,deltat)
% toc