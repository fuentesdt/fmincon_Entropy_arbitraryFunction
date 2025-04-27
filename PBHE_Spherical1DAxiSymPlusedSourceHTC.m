tic
clear; close all; clc;
%%% Solves the Pennes BioHeat Equation in Spherical Coordinates %%%
%%% The model is 1-Dimensional and Axisymmetric %%%
%%% The initial conditions are T(0,r) = 37 °C %%%
%%% Boundary Conditions are Symmetry at r = 0; and Fixed Temperature at r = R as T(0,R) = 37 °C %%%
%%% Uses TDMA algorithm to solve Discretized Linear Equations %%%
%%% Heat source by MNP is Pulsating in Time %%%
%% Problem Parameters
% debug
k_t=.4947; % [W/m/K] Thermal Conductivity
w_t = 0.0048; % [1/s] Blood perfusion rate
rhocp_t = 3557545; % [J/m^3/K] Density and Specific Heat Production
Qm_t = 0;                   % [W/m3] Tumor Metabolic Heat
% Md = 4.46e5; % [A/m] Domain Magnetization
Md = 7.77e5; % [A/m] Domain Magnetization

%%% Domain Parameters
RT = 0.005;               % [m] Radius of the Tumor
R = 0.0005+RT;                 % [m] Radius of the computaional domain
vol_t = (4*pi*RT^3)/3;    % [m3] Tumor Volume
dr = RT/50;               % [m] Spatial discretization step size
r = (0:dr:R)';            % [m] Spatial loactions
N = length(r);            % [-] Number of spatial nodes
t_loc = r<=RT;            % [-] Tumor r indices

%%% Time discretization: Implicit method is employed
t_end = 180; % [s] End time of the simulation
deltat = 3; % [s] Delta
Ntime = t_end/deltat; % [-] Number of time step
dt = 1; % Simulations time-step
t = 0:dt:t_end;
TS = length(t);

%%% Magnetic Field Parameters
Hmin = 7957;
Hmax = 39788;
% Hmin = 10000;
% Hmax = 10000;
H_range = Hmin:Hmin:Hmax;

% Generate random indices
indices = randi(length(H_range), 1, Ntime);

% Create the random array
Htime = H_range(indices);
for i = 1: length(Htime)
    % Calculate the starting index for the current iteration
    startIndex = (i - 1) * (deltat)/dt + 1;
    % Calculate the ending index for the current iteration
    endIndex = min(startIndex + (deltat)/dt - 1, TS-1);
    HVector(startIndex:endIndex) = Htime(i);
end

%%% Ambient Conditions
htc = 10;
T_amb = 24;

%%% Blood Properties
rho_b = 1000;           % [kg/m3] Blood Density
cp_b = 3840;            % [J/kg/K] Blood Specific Heat
k_b = 0.50;             % [W/m/K] Blood Thermal Conductivity
T_b = 37;               % [°C] Blood Arterial Temperature
T_bK = T_b + 273.15;    % [K] Blood Arterial Temperature

T_initial = 37;     % [°C] Initial Condition Temperature
alpha_t = k_t/(rhocp_t); % Tissue Thermal diffusivity [m^2/s]
f = 1e5;    % [Hz] Magnetic Field Frequency f = 100 kHz

%%% Magnetic Nanoparticle Physical and Magnetic Parameters MNP Used: Fe3O4
d_mnp = 6e-8;                 % [m] MNP Diameter d = 60 nm
delta_l = 1e-9;                 % [m] Liquid Layer Thickness on the Hard Solid MNP delta = 1 nm
dh_mnp = d_mnp+2*delta_l;       % [m] MNP Hydrodynamic Diameter
rho_mnp = 5180;                 % [kg/m3] MNP Density
Keff = 15000;                   % [J/m3] Effective Anisotropy Constant Keff = 20 kJ/m3

%%% MNP Dose and Magnetic Fluid Parameters
MNP_conc_Tumor = 1;                 % [mg/ml] MNP conc in injected fluid
MF_inj_vol = 0.1;                   % [ml, milliliters] Magnetic Fluid Volume injected into the Tumor
MF_vol = MF_inj_vol*1e-6;           % [m3] Magnetic Fluid Volume Injected in to the Tumor % Converion Factor for ml to m3 is 1e-6
mu_cf = 1e-3;                       % [Pa.s] Viscocity of the Carrier Fluid
MNP_mass = MNP_conc_Tumor*MF_inj_vol/1000;  % [kg] MNP mass in the fluid (=injected into the tumor); Varies with MNP Conc in MF
MNP_vol = MNP_mass/rho_mnp;         % [m3] volume of MNP in the Magnetic Fluid
MNP_vol_frac = MNP_vol/MF_vol;      % [-] Volume Fraction of MNP in MF

%%% Constants used in SAR Calculation
mu0 = pi*4e-7;              % [H/m] Permeability of Free Space
kB = 1.38e-23;              % [J/K] Blotzmann Constant
MNP_svol = pi*(d_mnp^3)/6;  % [m3] Solid Volume of Individual Mangetic Nanoparticle
MNP_hvol = pi*(dh_mnp^3)/6; % [m3] Hydrodynamic Volume of Individual Mangetic Nanoparticle
tau0 = 1e-9;                % [s] Attempt Time for Mangetic Moment relaxation
omega = 2*pi*f;             % [rad/s] Angular frequency of applied magnetic field
alpha_CF = 0.55;            % [-] Correction Factor for MNP Heating in MF and Tissue Mediums

%%% SAR Calculations
gamma = Keff*MNP_svol/(kB*T_bK);                    % [-] An intermediate parameters used in further calculations
tauB = (3*mu_cf*MNP_hvol)/(kB*T_bK);                % [s] Brownian Relaxation Time of MNP
tauN = sqrt(pi)*tau0*exp(gamma)/(2*sqrt(gamma));    % [s] Neel Relaxation Time of MNP
tauE = tauB*tauN/(tauB+tauN);                       % [s] Effective Relaxation Time of MNP

%% Implicit Finite Difference Scheme using the Cauchy Boundary Condition on the outer boundary
T = zeros(N, TS);
T(:,1) = T_initial; % Initial Condition
%%% Coefficient Matrix Three-Diagonal Generation
Lower = zeros(N, 1);
Main  = zeros(N, 1);
Upper = zeros(N, 1);
for i = 1:N
    if i ==1 % Symmetry Node
        Main(i) = 1 + 2*alpha_t*dt/(dr^2) + rho_b*cp_b*w_t*dt/(rhocp_t);
        Upper(i) = -2*alpha_t*dt/(dr^2);
    elseif i == N % Fixed Temperature Node
        Lower(i) = -k_t/dr;
        Main(i) = k_t/dr+htc;
    else % Internal Nodes
        Lower(i) = -alpha_t*dt/(dr^2) + alpha_t*dt/(r(i)*dr);
        Main(i) = 1 + 2*alpha_t*dt/(dr^2) + rho_b*cp_b*w_t*dt/(rhocp_t);
        Upper(i) = -alpha_t*dt/(dr^2) - alpha_t*dt/(r(i)*dr);
    end
end
%%% Solver Time Iterations
for n = 2:TS
    H = HVector(n-1); %
    if H == 0
        SAR = 0;
    else
        zeta = mu0*Md*H*MNP_svol/(kB*T_bK);                 % [-] Langevin parameter
        X_i = mu0*(Md^2)*MNP_vol_frac*MNP_svol/(3*kB*T_bK); % [-] Initial Magnetic Susceptibility
        X_0 = 3*X_i*(coth(zeta)-1/zeta)/zeta;               % [-] Equilibrium Magnetic Susceptibility
        X_L = X_0*omega*tauE/(1+(omega*tauE)^2);            % [-] Loss Component of Magnetic Susceptibility
        P = pi*mu0*f*X_L*H^2;                               % [W/m3] Heat Generation rate of MNPs in MF
        SAR = P/(rho_mnp*MNP_vol_frac);                     % [W/kg] MNP Specific Absorption Rate
        SAR_grams = SAR/1000;                               % [W/g] MNP Specific Absorption Rate
    end
    %%% Heat Source by MNP in Tumor
    MNP_conc_Tumor = MNP_mass/vol_t;      % [kg/m3] Concentration of MNP in Tumor, Assuming MNPs are distributed uniformly and confined within the tumor only
    Q_MNP = alpha_CF*MNP_conc_Tumor*SAR;  % [W/m3] Heat Generation by MNP in Tissue

    %%% Assign Q_MNP to the center region only
    q_mnp = Q_MNP*t_loc;        % MNP heat generation [W/m^3] On/Off Pulsating with time
    q_mnp_time(:,n-1) = q_mnp;

    Force = zeros(N,1); % Force vector
    Force(1:N-1)  = T(1:N-1,n-1) + (rho_b*cp_b.*w_t*dt*T_b + Qm_t*dt+ q_mnp(1:N-1)*dt)/(rhocp_t);
    Force(N) = htc*T_amb;
    T(:,n) = thomas_algorithm(Lower, Main, Upper, Force); % Spatio-Temporal Temperature
end

%%
T_sum_time = sum(T,2);      % Summation of Temperature over all Times at each r location
G = trapz(r, T_sum_time);     % Integration of Temperature over the domain in 3D
tempqoi = G;
%% Plot the results
q_mnp_t = [q_mnp_time(:,1),q_mnp_time];
[RR,TT] = meshgrid(r,t);
figure('Position', [40, 40, 800, 600])
subplot(2,2,1)
plot((1:length(HVector))*dt,HVector)
xlim([0,t(end)])
xticks([0:t(end)/4:t(end)]);
xlabel('Time [s]')
ylabel('H Amplitude, [A/m]');
title('Magnetic Field');

subplot(2,2,2)
surf(RR,TT,q_mnp_t')
xlabel('Radial Distance [m]')
ylabel('Time [s]');
zlabel('Qmnp [W/m^3]');
title('Spatio-Temporal Heat Source');
xlim([0,R])
xticks([0,RT/2,RT,R]);
ylim([0,t(end)])
yticks([0:t(end)/4:t(end)]);
grid on;
subplot(2,2,3)
plot_radius = [0,RT,R];
legend_String = string(plot_radius)+[" m Tumor Center"," m Tumor Edge"," m Outer Boundary"];
plot_rad_idx = floor((plot_radius/dr))+1;
plot(t,T(plot_rad_idx,:), 'LineWidth',2)
xlabel('Time, t [s]');
ylabel('Temperature, T [°C]');
xlim([0,t(end)])
xticks([0:t(end)/4:t(end)]);
legend(legend_String, Location='best')
title('Temperature Elevations');
grid on;
subplot(2,2,4)
plot_time = [60,120,t(end)];
legend_String = string(plot_time)+[" s"," s"," s"];
plot_ind = (plot_time/dt)+1;
plot(r,T(:,plot_ind),'LineWidth',2)
hold on
xlabel('Radial distance, r [m]');
xlim([0,R])
xticks(0:R/5:R);
ylabel('Temperature, T [°C]');
legend(legend_String, Location='best')
title('Spatial Temperature profile');
grid on;

%% Tri-Diagonal Matrix Alogorithm for Linear System of Equation Solver
function x = thomas_algorithm(Lower, Main, Upper, Force)
nn = length(Force);
c_star = zeros(nn, 1);
d_star = zeros(nn, 1);
x = zeros(nn, 1);

c_star(1) = Upper(1) / Main(1);
d_star(1) = Force(1) / Main(1);

for ii = 2:nn
    c_star(ii) = Upper(ii) / (Main(ii) - Lower(ii) * c_star(ii-1));
    d_star(ii) = (Force(ii) - Lower(ii) * d_star(ii-1)) / (Main(ii) - Lower(ii) * c_star(ii-1));
end

x(nn) = d_star(nn);
for ii = nn-1:-1:1
    x(ii) = d_star(ii) - c_star(ii) * x(ii+1);
end
end