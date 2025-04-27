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

    for j = 1:size(data,1)
        Qm_t = 0;                   % [W/m3] Tumor Metabolic Heat
        k_t=data(j,1); % [W/m/K] Thermal Conductivity
        w_t = data(j,2); % [1/s] Blood perfusion rate
        rhocp_t = data(j,3); % [J/m^3/K] Density and Specific Heat Production
        Km = data(j,4); % [J/m^3] Anisotropy Constant
        Gain(j,i) = pennesmht_fe(k_t,w_t,rhocp_t,Km,Htime,deltat);
    end
end
%%
close all
figure()
for i = 1:size(Gain,2)
subplot(2,4,i)
plot(Gain(:,i))
xlabel('Qudrature Points')
ylabel('Gain G(K,P) @ fixed K')
title(['Initial Guess [K] :',num2str(i)])
end
%%

for i = 1:size(Gain,1)
    figure()
plot(Gain(i,:))
%pause
xlabel('Initiall Guess [K]')
ylabel('Gain G(K,P)')
title(['QP [P] :',num2str(i)])
saveas(gca,['Gainforqp',num2str(i)],'png')
close
end

%%
figure()
histogram(Gain)

xlabel('Gain Value')
ylabel('Frequncey')

%%
meangain = mean(Gain,1);
figure()
plot(meangain)
xlabel('Initial Guesses')
ylabel('Mean G(K,P) over 81 samples')

%%
close all
figure()
for i = 1:size(Gain,2)
subplot(2,4,i)
histogram(Gain(:,i))
% xlabel('Qudrature Points')
% ylabel('Gain G(K,P) @ fixed K')
% title(['Initial Guess [K] :',num2str(i)])
end
