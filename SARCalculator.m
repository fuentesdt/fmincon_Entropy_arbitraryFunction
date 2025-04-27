clc
clear
close all
T_medium = 37;               % [°C] Blood Arterial Temperature
T_mediaumK = T_medium + 273.15;    % [K] Blood Arterial Temperature

f = 1.60e5;    % [Hz] Magnetic Field Frequency f = 163 kHz
H = 7600;
%%% Magnetic Nanoparticle Physical and Magnetic Parameters MNP Used: Fe3O4
d_mnp = 70e-9;                 % [m] MNP Diameter d = 60 nm
delta_l = 0.1*d_mnp;                 % [m] Liquid Layer Thickness on the Hard Solid MNP delta = 1 nm
dh_mnp = d_mnp+2*delta_l;       % [m] MNP Hydrodynamic Diameter

rho_mnp = 5180;        % [kg/m3] MNP Density
% Md_measured = 150;                  % [emu/g]
% Md = Md_measured*rho_mnp;         % [A/m] Domain Magnetization
Md = 446000;         % [A/m] Domain Magnetization
Km = 23;
%%% MNP Dose and Magnetic Fluid Parameters
mnp_mf_conc = 1000; % [ug/ml] micrograms of MNP/milliliter of water
mnp_mf_conc2 = mnp_mf_conc/1000; %[kg/m^3] kg of MNP/m3 of water
mf_vol = 100; %[ul] microliters of nanofluid injected
mf_vol2 = mf_vol*1e-9; % [m3] of nanofluid
mnp_mass_kg = mnp_mf_conc2*mf_vol2;
mnp_vol = mnp_mass_kg/(rho_mnp); % [m3] of MNP
mu_cf = 1e-3;                   % [Pa.s] Viscocity of the Carrier Fluid
mnp_vol_frac = mnp_vol/mf_vol2;  %[-] Volume Fraction of MNP in MF

%%% Constants used in SAR Calculation
mu0 = pi*4e-7;              % [H/m] Permeability of Free Space
kB = 1.38e-23;              % [J/K] Blotzmann Constant
mnp_svol = pi*(d_mnp^3)/6;  % [m3] Solid Volume of Individual Mangetic Nanoparticle
mnp_hvol = pi*(dh_mnp^3)/6; % [m3] Hydrodynamic Volume of Individual Mangetic Nanoparticle
tau0 = 1e-9;                % [s] Attempt Time for Mangetic Moment relaxation
omega = 2*pi*f;             % [rad/s] Angular frequency of applied magnetic field
alpha_CF = 0.55;            % [-] Correction Factor for MNP Heating in MF and Tissue Mediums

%%% SAR Calculations
gamma = Km*mnp_svol/(kB*T_mediaumK);                    % [-] An intermediate parameters used in further calculations
tauB = (3*mu_cf*mnp_hvol)/(kB*T_mediaumK);                % [s] Brownian Relaxation Time of MNP
tauN = tau0*exp(gamma);    % [s] Neel Relaxation Time of MNP
% if isinf(tauN)
%     error(['tauN is infinity and gamma = ', num2str(gamma)])
% end
One_by_tauE = 1/tauB+1/tauN;                       % [s] Effective Relaxation Time of MNP
tauE = 1/One_by_tauE;

zeta = mu0*Md*H*mnp_svol/(kB*T_mediaumK);                 % [-] Langevin parameter
X_i = mu0*(Md^2)*mnp_vol_frac*mnp_svol/(3*kB*T_mediaumK); % [-] Initial Magnetic Susceptibility
X_0 = 3*X_i*(coth(zeta)-1/zeta)/zeta;               % [-] Equilibrium Magnetic Susceptibility
X_L = X_0*omega*tauE/(1+(omega*tauE)^2);            % [-] Loss Component of Magnetic Susceptibility
P = pi*mu0*f*X_L*H^2;                               % [W/m3] Heat Generation rate of MNPs in MF
SAR = P/(rho_mnp*mnp_vol_frac);                     % [W/kg] MNP Specific Absorption Rate
SAR_grams = SAR/1000                              % [W/g] MNP Specific Absorption Rate
