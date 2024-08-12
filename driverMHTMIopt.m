clear all
clc

%myoptions.Algorithm = 'constDirect'
%myoptions = optimoptions(@fmincon,'Display','iter-detailed','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'MaxFunctionEvaluations',1e7,'ConstraintTolerance',2.e-9, 'OptimalityTolerance',2.5e-4,'Algorithm','interior-point','StepTolerance',1.000000e-12,'MaxIterations',1000,'PlotFcn',{'optimplotfvalconstr', 'optimplotconstrviolation', 'optimplotfirstorderopt' },'HonorBounds',true, 'Diagnostic','on','FunValCheck','on' )
myoptions = optimoptions(@fmincon,'Display','iter-detailed','MaxFunctionEvaluations',1e7,'ConstraintTolerance',2.e-9, 'OptimalityTolerance',2.5e-4,'Algorithm','sqp','StepTolerance',1.000000e-12,'MaxIterations',1000,'PlotFcn',{'optimplotfvalconstr', 'optimplotconstrviolation', 'optimplotfirstorderopt' },'HonorBounds',true, 'Diagnostic','on','FunValCheck','on','FiniteDifferenceStepSize',1 )

%% % monitor memory: while [ -e /proc/3291925 ] ; do  top -b -n 1 -p 3291925 >>process.txt ;sleep 60; done

%function driverHPMIopt(NGauss,NumberUncertain,modelSNR,myoptions,ObjectiveType,GaussLegendre )
NGauss = 3,NumberUncertain=1,modelSNR=10, ObjectiveType = 'TotalSignal',GaussLegendre=false

NGauss,NumberUncertain,modelSNR,myoptions.Algorithm,ObjectiveType,GaussLegendre
close all


%% Variable Setup
t_end = 180; % [s] End time of the simulation
deltat = 3; % [s] Delta
Ntime = t_end/deltat; % [-] Number of time step
Nspecies = 1
% switch between uniform and Gaussian RV for prior
if(GaussLegendre)
    QuadratureRule = 'Legendre';
else
    QuadratureRule = 'Hermite';
end

% debug
k_t=.527; % [W/m/K] Thermal Conductivity
w_t = 0.009; % [1/s] Blood perfusion rate
rhocp_t = 1045*3600; % [J/m^3/K] Density and Specific Heat Production
Md = 4.46e5; % Domain Magnetization

%%% Time varying Magnetic Field Amplitude Vector Creation
Hmin = 7957; % [A/m] Minimum value of Magnetic Field the device can generate
Hmax = 39788; % [A/m] Maximum value of Magnetic Field the device can generate
H_range = (Hmin:Hmin:Hmax)'; % Possible steps of Magnetic Field on the device
indices = randi(length(H_range), 1, Ntime); % Generate random indices

Htime = H_range(indices); % Random variation of H in time and constant within each step of "deltat"
Gain_Variable_H = pennesmht(k_t,w_t,rhocp_t,Md,Htime,deltat);


Htime = 10000*ones(Ntime,1); % Constant H in time
tic
Gain_H_10000 = MIGHQuadMHT(Htime,NGauss,NumberUncertain,Nspecies,Ntime,GaussLegendre,ObjectiveType,deltat);
toc

FDHtime = Htime;
deltaH = 100;
FDHtime(10) = FDHtime(10)+deltaH; 

Gain_H_10000FD = MIGHQuadMHT(FDHtime,NGauss,NumberUncertain,Nspecies,Ntime,GaussLegendre,ObjectiveType,deltat);

% @mahesh why is this FD gradient so small ? can you debug ? 
(Gain_H_10000FD  -Gain_H_10000  )
(Gain_H_10000FD  -Gain_H_10000  )/deltaH

%% optimize MI
optf = false;
if optf

    tic;

    %%
    % Create an optimization problem using these converted optimization expressions.
    disp('create optim prob')

    % compare walker solution at qp
    switch (NumberUncertain)
        case(3)
            optparams = struct('t0',[t0qp(1);0],'gammaPdfA',[alphamean(1)  ;1],'gammaPdfB',[betamean(1);1],...
                'scaleFactor',VIF_scale_fact,'T1s',[T1pmean(1),T1lmean(1)],'ExchangeTerms',[0,kplqp(1) ;0,0],...
                'TRList',TR_list,'PerfusionTerms',[kveqp(1),0],'volumeFractions',ve,...
                'fitOptions', opts);
        case(4)
            error("WIP")
        case(5)
            error("WIP")
    end
    % solve
    switch (myoptions.Algorithm)
        case('constDirect')
            powergrid = [0:1000:40000];
            brutesearch= zeros(size(powergrid ));
            backspaces = '';
            for iii = 1:length(powergrid(:))
                %disp(sprintf('%d %d ',pyrgrid(iii),lacgrid(iii)));
                % Print percentage progress
                percentage = iii/length(powergrid(:));
                perc_str = sprintf('completed %3.1f', percentage);
                fprintf([backspaces, perc_str]);
                backspaces = repmat(sprintf('\b'), 1, length(perc_str));
                constH = powergrid(iii)*  ones(Ntime,1);
                brutesearch(iii) = MIGHQuadMHT(constH,NGauss,NumberUncertain,Nspecies,Ntime,GaussLegendre,ObjectiveType,deltat);
            end
            save(sprintf('brutesearchNG%dNu%d%s%sSNR%02d%s.mat',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule) ,'brutesearch','powergrid')
            [maxMI,idmax] = max(brutesearch(:));
            [fval,idmin] = min(brutesearch(:));
            popt.powerList = repmat(powergrid(idmin),1,Ntime);
            %pneg.FaList = repmat([pi/180*pyrgrid(idmax);pi/180*lacgrid(idmax)],1,Ntime);
            %optparams.FaList = repmat([pi/180*pyrgrid(idmin);pi/180*lacgrid(idmin)],1,Ntime);
            %refparams.FaList = repmat([pi/180*pyrgrid(idmin);pi/180*lacgrid(idmin)],1,Ntime);

            handle = figure(5)
            plot(powergrid,brutesearch)
            xlabel('power ')
            ylabel('MI ')
            title(sprintf('max %f min %f',maxMI,fval) )
            %text(pyrgrid(idmin)+1,lacgrid(idmin)+1, sprintf('opt %d %d', pyrgrid(idmin), lacgrid(idmin)));
            %text(pyrgrid(idmax)+1,lacgrid(idmax)+1, 'control');
        otherwise
            InitialGuess = 10000*rand(Ntime,1); 
            pmin =  [Htime(:)*0];
            pmax =  [Htime(:)*0+40000];
            tolx=1.e-9;
            tolfun=5.e-4;
            maxiter=400;

            % truthconstraint = infeasibility(stateconstraint,x0);
            %[popt,fval,exitflag,output] = solve(convprob,x0,'Options',myoptions, 'ConstraintDerivative', 'auto-reverse', 'ObjectiveDerivative', 'auto-reverse' )
            Fx = @(x) MIGHQuadMHT(x,NGauss,NumberUncertain,Nspecies,Ntime,GaussLegendre,ObjectiveType,deltat);
            [designopt,fval,exitflag,output,lambda,grad,hessian] ...
                =fmincon(Fx, InitialGuess ,[],[],[],[],pmin,pmax,[],myoptions);

            handle = figure(5)
            %% optparams.FaList = reshape(designopt(:),size(params.FaList ));
            %% refparams.FaList = reshape(designopt(:),size(params.FaList ));
            %% popt.FaList      = reshape(designopt(:),size(params.FaList ));
    end
    % save convergence history
    set(gca,'FontSize',16)
    saveas(handle,sprintf('historyNG%dNu%d%s%sSNR%02d%s',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule ),'png')
    toc;


    %% [t_axisopt,Mxyopt,Mzopt] = model.compile(M0.',optparams);
    %% [t_axisref,Mxyref,Mzref] = model.compile(M0.',refparams);
    %% save(sprintf('poptNG%dNu%d%s%sSNR%02d%s.mat',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule) ,'fval','popt','params','Mxy','Mz','Mxyref','Mzref','signu','signuImage')
    %% handle = figure(10)
    %% plot(params.TRList,Mxyopt(1,:),'b',params.TRList,Mxyopt(2,:),'k')
    %% ylabel('MI Mxy')
    %% xlabel('sec'); legend('Pyr','Lac')
    %% set(gca,'FontSize',16)
    %% saveas(handle,sprintf('OptMxyNG%dNu%d%s%sSNR%02d%s',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule),'png')
    %% handle = figure(11)
    %% plot(params.TRList,popt.FaList(1,:)*180/pi,'b',params.TRList,popt.FaList(2,:)*180/pi,'k')
    %% ylabel('MI FA (deg)')
    %% xlabel('sec'); legend('Pyr','Lac')
    %% ylim([0 40])
    %% set(gca,'FontSize',16)
    %% saveas(handle,sprintf('OptFANG%dNu%d%s%sSNR%02d%s',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule),'png')
    %% handle = figure(12)
    %% plot(params.TRList,Mzopt(1,:),'b--',params.TRList,Mzopt(2,:),'k--')
    %% hold
    %% plot(params.TRList,cos(optparams.FaList(1,:))'.*popt.state(:,1, 1),'b',params.TRList,cos(optparams.FaList(2,:))'.*popt.state(:,2, 1),'k')
    %% %if(lqp > 1)
    %% %  plot(params.TRList,popt.state(:,1, 5),'b',params.TRList,popt.state(:,2, 5),'k')
    %% %  plot(params.TRList,popt.state(:,1,10),'b',params.TRList,popt.state(:,2,10),'k')
    %% %  plot(params.TRList,popt.state(:,1,15),'b',params.TRList,popt.state(:,2,15),'k')
    %% %end
    %% ylabel('MI Mz ')
    %% xlabel('sec'); legend('Pyr','Lac')
    %% set(gca,'FontSize',16)
    %% saveas(handle,sprintf('OptMzNG%dNu%d%s%sSNR%02d%s',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule),'png')
end


% evaluate MI
function MIobjfun =MIGHQuadMHT(hOpt,NGauss,NumberUncertain,Nspecies,Ntime,GaussLegendre,ObjectiveType,deltat)
    %% Thermal conductivity 
    kmean = [ .5 ]; % s
    kstdd = [ .2 ]; % s
    klb   = [ 0.1  ]; % s
    kub   = [ 0.8 ]; % s
    %% perfusion 
    wmean = [ 0.006  ]; % s
    wstdd = [ 0.003  ]; % s
    wlb   = [ 0.002  ]; % s
    wub   = [ 0.010 ]; % s
    %% rho * specific heat 
    crhomean = [ 3762000 ];       % s
    crhostdd = [ 100000 ];       % s
    crholb   = [ 3600000 ];       % s
    crhoub   = [ 4000000 ];       % s
    %% domain magnetization 
    mdmean = [ 446000 ];      % s
    mdstdd = [ 100000  ];      % s
    mdlb   = [ 200000 ];      % s
    mdub   = [ 500000 ];      % s

    %% signal uncertianty
    signu = sqrt(2* Ntime) * .1;
    [x2,xn2,xm2,w2,wn2]=GaussHermiteNDGauss(NGauss,0,signu);
    lqp2=length(xn2{1}(:));

    % switch between uniform and Gaussian RV for prior
    if(GaussLegendre)
        QuadratureRule = 'Legendre';
    else
        QuadratureRule = 'Hermite';
    end
    %disp('evaluate quadrature points')
    switch (NumberUncertain)
        case(1)
            if(GaussLegendre)
                [x,xn,xm,w,wn]=GaussLegendreNDGauss(NGauss,mdlb,mdub);
            else
                [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,mdmean,mdstdd);
            end
            mdqp   = xn{1}(:);
            % evaluate function at each quadrature point
            lqp=length(xn{1}(:));
            sumstatevariable      = zeros(Nspecies,lqp);
            for iqp = 1:lqp
              sumstatevariable(1,iqp) =  pennesmht(kmean,wmean,crhomean,mdqp(iqp),hOpt,deltat);
            end
        case(4)
            if(GaussLegendre)
                [x,xn,xm,w,wn]=GaussLegendreNDGauss(NGauss,[klb; wlb; crholb; mdlb],[kub; wub; crhoub; mdub]);
            else
                [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,[kmean; wmean; crhomean; mdmean],[kstdd; wstdd; crhostdd; mdstdd]);
            end
            kqp    = xn{1}(:);
            wqp    = xn{2}(:);
            crhoqp = xn{3}(:);
            mdqp  = xn{4}(:);
            % evaluate function at each quadrature point
            lqp=length(xn{1}(:));
            sumstatevariable      = zeros(Nspecies,lqp);
            for iqp = 1:lqp
              sumstatevariable(1,iqp) =  pennesmht(kqp(iqp),wqp(iqp),crhoqp(iqp),mdqp(iqp),hOpt,deltat);
            end
    end

    %disp('build objective function')
    expandvar  = ones(1,lqp);

    switch (ObjectiveType)
        case('TotalSignal')
            diffsumm =sumstatevariable(1,:)' * expandvar   - expandvar' * sumstatevariable(1,:);
            negHz = 0;
            for jjj=1:lqp2
                znu=xn2{1}(jjj) ;
                % note that the sqrt(pi)^N+1  and 2^N factors from the integration over priors is included in the quadrature routines.
                % this makes is easier to switch between Gaussian and Uniform RV
                negHz = negHz + wn2(jjj) * (wn(:)' * log(exp(-(znu + diffsumm).^2/2/signu^2 - log(signu) -log(2*pi)/2   ) * wn(:)));
            end
    end
    % MI = H(z) - H(z|P)
    %  H(z|P)  constant ==> max MI = max H(z) = min -H(z)
    MIobjfun = negHz;
end

function tempqoi = pennesmht(k_t,w_t,rhocp_t,Md,Htime,deltat)
%%% Following are the inputs for the function
% rhocp_t = Tissue Density [kg/m3] * cp_t = Tissue Specific Heat [J/kg/K]
% k_t = Tissue Thermal Conductivity [W/m/K]
% w_t = Blood Perfusion Rate [1/s]
% Qm_t = Metabolic Heat Generation Rate [W/m3]
% H = Applied Magnetic Field [A/m]
% Md = Domain Magnetization Value of MNP [A/m]
% Keff = Magnetic Anistropy Constant of MNP [J/m3]
% deltat = time step in which the H is constant [s]
Qm_t = 0;                   % [W/m3] Tumor Metabolic Heat

%%% Domain Parameters
RT = 0.005;               % [m] Radius of the Tumor
R = 1.1*RT;               % [m] Radius of the computaional domain
vol_t = (4*pi*RT^3)/3;    % [m3] Tumor Volume
dr = RT/50;               % [m] Spatial discretization step size
r = (0:dr:R)';            % [m] Spatial loactions
N = length(r);            % [-] Number of spatial nodes
t_loc = r<=RT;            % [-] Tumor r indices

%%% Time discretization: Implicit method is employed
dt = 1; % time-step for forward runs
t = 0:size(Htime,1)*deltat;     % Discrete times
TS = length(t);     % Number of time steps

%%% Assigning the PieceWise Stepped H in deltat at respective times at dt resolution
for i = 1: length(Htime)
    startIndex = (i - 1) * (deltat)/dt + 1; % Calculate the starting index for the current iteration
    endIndex = min(startIndex + (deltat)/dt - 1, TS-1); % Calculate the ending index for the current iteration
    HVector(startIndex:endIndex) = Htime(i); % Amplitude vector in time at a resolution of dt assigned from deltat resolution  
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
d_mnp = 1.6e-8;                 % [m] MNP Diameter d = 16 nm
delta_l = 1e-9;                 % [m] Liquid Layer Thickness on the Hard Solid MNP delta = 1 nm
dh_mnp = d_mnp+2*delta_l;       % [m] MNP Hydrodynamic Diameter
rho_mnp = 5180;                 % [kg/m3] MNP Density
Keff = 20000;                   % [J/m3] Effective Anisotropy Constant Keff = 20 kJ/m3

%%% MNP Dose and Magnetic Fluid Parameters
MNP_dose = 4;                   % [kg/m3] MNP Dose MNP mass used per unit vol of tumor = 4 mg/cm3; % Conversion 1 mg/cm3 = 1 kg/m3
MNP_mass = MNP_dose*vol_t;      % [kg] MNP mass injected into the tumor; Varies with MNP Dose
MNP_vol = MNP_mass/rho_mnp;     % [m3] volume of MNP in the Magnetic Fluid
MF_inj_vol = 1;                 % [ml, milliliters] Magnetic Fluid Volume injected into the Tumor
MF_vol = MF_inj_vol*1e-6;       % [m3] Magnetic Fluid Volume Injected in to the Tumor % Converion Factor for ml to m3 is 1e-6
mu_cf = 1e-3;                   % [Pa.s] Viscocity of the Carrier Fluid
MNP_vol_frac = MNP_vol/MF_vol;  %[-] Volume Fraction of MNP in MF

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
    zeta = mu0*Md*H*MNP_svol/(kB*T_bK);                 % [-] Langevin parameter
    X_i = mu0*(Md^2)*MNP_vol_frac*MNP_svol/(3*kB*T_bK); % [-] Initial Magnetic Susceptibility
    X_0 = 3*X_i*(coth(zeta)-1/zeta)/zeta;               % [-] Equilibrium Magnetic Susceptibility
    X_L = X_0*omega*tauE/(1+(omega*tauE)^2);            % [-] Loss Component of Magnetic Susceptibility
    P = pi*mu0*f*X_L*H^2;                               % [W/m3] Heat Generation rate of MNPs in MF
    SAR = P/(rho_mnp*MNP_vol_frac);                     % [W/kg] MNP Specific Absorption Rate
    SAR_grams = SAR/1000;                               % [W/g] MNP Specific Absorption Rate

    %%% Heat Source by MNP in Tumor
    MNP_conc = MNP_mass/vol_t;      % [kg/m3] Concentration of MNP in Tumor, Assuming MNPs are distributed uniformly and confined within the tumor only
    Q_MNP = alpha_CF*MNP_conc*SAR;  % [W/m3] Heat Generation by MNP in Tissue

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
%% j%% Plot the results
%% jq_mnp_t = [q_mnp_time(:,1),q_mnp_time];
%% j[RR,TT] = meshgrid(r,t);
%% jfigure('Position', [40, 40, 800, 600])
%% jsubplot(2,2,1)
%% jplot((1:length(HVector))*dt,HVector)
%% jxlim([0,t(end)])
%% jxticks([0:t(end)/4:t(end)]);
%% jxlabel('Time [s]')
%% jylabel('H Amplitude, [A/m]');
%% jtitle('Magnetic Field');
%% j
%% jsubplot(2,2,2)
%% jsurf(RR,TT,q_mnp_t')
%% jxlabel('Radial Distance [m]')
%% jylabel('Time [s]');
%% jzlabel('Qmnp [W/m^3]');
%% jtitle('Spatio-Temporal Heat Source');
%% jxlim([0,R])
%% jxticks([0,RT/2,RT,R]);
%% jylim([0,t(end)])
%% jyticks([0:t(end)/4:t(end)]);
%% jgrid on;
%% j
%% jsubplot(2,2,3)
%% jplot_radius = [0,RT,R];
%% jlegend_String = string(plot_radius)+[" m Tumor Center"," m Tumor Edge"," m Outer Boundary"];
%% jplot_rad_idx = (plot_radius/dr)+1;
%% jplot(t,T(plot_rad_idx,:), 'LineWidth',2)
%% jxlabel('Time, t [s]');
%% jylabel('Temperature, T [°C]');
%% jxlim([0,t(end)])
%% jxticks([0:t(end)/4:t(end)]);
%% jlegend(legend_String, Location='best')
%% jtitle('Temperature Elevations');
%% jgrid on;
%% j
%% jsubplot(2,2,4)
%% jplot_time = [60,120,t(end)];
%% jlegend_String = string(plot_time)+[" s"," s"," s"];
%% jplot_ind = (plot_time/dt)+1;
%% jplot(r,T(:,plot_ind),'LineWidth',2)
%% jhold on
%% jxlabel('Radial distance, r [m]');
%% jxlim([0,R])
%% jxticks(0:R/5:R);
%% jylabel('Temperature, T [°C]');
%% jlegend(legend_String, Location='best')
%% jtitle('Spatial Temperature profile');
%% jgrid on;

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
end
