clear
clc
set(0, 'DefaultAxesXColor', 'k', ...
       'DefaultAxesYColor', 'k', ...
       'DefaultTextColor', 'k');
% Define the number of workers (processors)
numWorkers = 10; % Replace 4 with the desired number of processors

% Start a parallel pool with the specified number of workers
if isempty(gcp('nocreate')) % Check if a parallel pool exists
    parpool(numWorkers); % Start a parallel pool with the specified number of workers
end

%myoptions.Algorithm = 'constDirect'
myoptions = optimoptions(@fmincon,'UseParallel', true, 'Display','iter-detailed','MaxFunctionEvaluations',1e7,'ConstraintTolerance',1e-20, 'OptimalityTolerance',1e-20,'Algorithm','sqp','StepTolerance',1e-20,'MaxIterations',400,'PlotFcn',{@optimplotfvalconstr, @optimplotconstrviolation, @optimplotfirstorderopt },'HonorBounds',true, 'Diagnostic','on','FunValCheck','on','FiniteDifferenceStepSize',1 )

%% % monitor memory: while [ -e /proc/3291925 ] ; do  top -b -n 1 -p 3291925 >>process.txt ;sleep 60; done

%function driverHPMIopt(NGauss,NumberUncertain,modelSNR,myoptions,ObjectiveType,GaussLegendre )
NGauss = 3,NumberUncertain=4,modelSNR=10, ObjectiveType = 'TotalSignal',GaussLegendre=false

NGauss,NumberUncertain,modelSNR,myoptions.Algorithm,ObjectiveType,GaussLegendre
close all


%% Variable Setup
t_end = 1800; % [s] End time of the simulation
deltat = 60; % [s] Delta
Ntime = t_end/deltat; % [-] Number of time step
Nspecies = 1
% switch between uniform and Gaussian RV for prior
if(GaussLegendre)
    QuadratureRule = 'Legendre';
else
    QuadratureRule = 'Hermite';
end

%%% Mean uncertain variables
k_t=.492; % [W/m/K] Thermal Conductivity
w_t = 0.00682; % [1/s] Blood perfusion rate
rhocp_t = 3589229; % [J/m^3/K] Density and Specific Heat Production
Km = 31500; % [J/m^3] Anisotropy Constant %% A variable to be recoved from the experiments/data

%%% Time varying Magnetic Field Amplitude Vector Creation
Hmin = 8000; % [A/m] Minimum value of Magnetic Field the device can generate
Hmax = 40000; % [A/m] Maximum value of Magnetic Field the device can generate
H_range = (Hmin:Hmin:Hmax)'; % Possible steps of Magnetic Field on the device

% %%% Calculating initial gain/total signal and the Mutual Information for all the trail initial guesses
% InitialGuess = 0*ones(Ntime,1); IG =1;
% Htime = InitialGuess; % Constant H in time
% eval(['Initial_Gain_IG',num2str(IG),' = pennesmht(k_t,w_t,rhocp_t,Km,Htime,deltat)']); % Intial gain = total signal for the IG
% eval(['MI_IG',num2str(IG),' = MIGHQuadMHT(Htime,NGauss,NumberUncertain,Nspecies,Ntime,GaussLegendre,ObjectiveType,deltat)']); %Initial MI for the IG
% 
% InitialGuess = H_range(1)*ones(Ntime,1);  IG =2;
% Htime = InitialGuess; % Constant H in time
% eval(['Initial_Gain_IG',num2str(IG),' = pennesmht(k_t,w_t,rhocp_t,Km,Htime,deltat)']); % Intial gain = total signal for the IG
% eval(['MI_IG',num2str(IG),' = MIGHQuadMHT(Htime,NGauss,NumberUncertain,Nspecies,Ntime,GaussLegendre,ObjectiveType,deltat)']); %Initial MI for the IG
% 
% InitialGuess = H_range(2)*ones(Ntime,1);  IG =3;
% Htime = InitialGuess; % Constant H in time
% eval(['Initial_Gain_IG',num2str(IG),' = pennesmht(k_t,w_t,rhocp_t,Km,Htime,deltat)']); % Intial gain = total signal for the IG
% eval(['MI_IG',num2str(IG),' = MIGHQuadMHT(Htime,NGauss,NumberUncertain,Nspecies,Ntime,GaussLegendre,ObjectiveType,deltat)']); %Initial MI for the IG
% 
% InitialGuess = H_range(3)*ones(Ntime,1);  IG =4;
% Htime = InitialGuess; % Constant H in time
% eval(['Initial_Gain_IG',num2str(IG),' = pennesmht(k_t,w_t,rhocp_t,Km,Htime,deltat)']); % Intial gain = total signal for the IG
% eval(['MI_IG',num2str(IG),' = MIGHQuadMHT(Htime,NGauss,NumberUncertain,Nspecies,Ntime,GaussLegendre,ObjectiveType,deltat)']); %Initial MI for the IG
% 
% InitialGuess = H_range(4)*ones(Ntime,1);  IG =5;
% Htime = InitialGuess; % Constant H in time
% eval(['Initial_Gain_IG',num2str(IG),' = pennesmht(k_t,w_t,rhocp_t,Km,Htime,deltat)']); % Intial gain = total signal for the IG
% eval(['MI_IG',num2str(IG),' = MIGHQuadMHT(Htime,NGauss,NumberUncertain,Nspecies,Ntime,GaussLegendre,ObjectiveType,deltat)']); %Initial MI for the IG
% 
% InitialGuess = H_range(5)*ones(Ntime,1);  IG =6;
% Htime = InitialGuess; % Constant H in time
% eval(['Initial_Gain_IG',num2str(IG),' = pennesmht(k_t,w_t,rhocp_t,Km,Htime,deltat)']); % Intial gain = total signal for the IG
% eval(['MI_IG',num2str(IG),' = MIGHQuadMHT(Htime,NGauss,NumberUncertain,Nspecies,Ntime,GaussLegendre,ObjectiveType,deltat)']); %Initial MI for the IG
% 
% indices = [3 5 1 2 4 2 3 5 1 3 2 4 5 2 1 4 5 3 1 5 4 2 3 1 5 2 3 1 5 2];
% InitialGuess = H_range(indices); IG =7; % 1. Random variation of H bounded by the minimum and maximum values
% Htime = InitialGuess; % Constant H in time
% eval(['Initial_Gain_IG',num2str(IG),' = pennesmht(k_t,w_t,rhocp_t,Km,Htime,deltat)']); % Intial gain = total signal for the IG
% eval(['MI_IG',num2str(IG),' = MIGHQuadMHT(Htime,NGauss,NumberUncertain,Nspecies,Ntime,GaussLegendre,ObjectiveType,deltat)']); %Initial MI for the IG
% 
% [H_var_steps,] = H_variation_between_H_AND_0(Ntime);
% InitialGuess = 7600*H_var_steps; IG =8;   % Rob's group H variation Peak amplitude was 7600 [A/m] with 50% duty cycle (60s On and 60s Off)
% Htime = InitialGuess; % Constant H in time
% eval(['Initial_Gain_IG',num2str(IG),' = pennesmht(k_t,w_t,rhocp_t,Km,Htime,deltat)']); % Intial gain = total signal for the IG
% eval(['MI_IG',num2str(IG),' = MIGHQuadMHT(Htime,NGauss,NumberUncertain,Nspecies,Ntime,GaussLegendre,ObjectiveType,deltat)']); %Initial MI for the IG


%% optimize MI
optf = true;
if optf

    tic;

    %%
    % Create an optimization problem using these converted optimization expressions.
    disp('create optim prob')

    % compare walker solution at qp
    % switch (NumberUncertain)
    %     case(3)
    %         optparams = struct('t0',[t0qp(1);0],'gammaPdfA',[alphamean(1)  ;1],'gammaPdfB',[betamean(1);1],...
    %             'scaleFactor',VIF_scale_fact,'T1s',[T1pmean(1),T1lmean(1)],'ExchangeTerms',[0,kplqp(1) ;0,0],...
    %             'TRList',TR_list,'PerfusionTerms',[kveqp(1),0],'volumeFractions',ve,...
    %             'fitOptions', opts);
    %     case(4)
    %         error("WIP")
    %     case(5)
    %         error("WIP")
    % end
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
            %%% These are traial initial guesses and should be used one at a time to find Optimum H
            clear InitialGuess IG
            % IG =1; InitialGuess = 0*ones(Ntime,1); 
            % IG =2; InitialGuess = H_range(1)*ones(Ntime,1);  
            % IG =3; InitialGuess = H_range(2)*ones(Ntime,1);  
            % IG =4; InitialGuess = H_range(3)*ones(Ntime,1);  
            % IG =5; InitialGuess = H_range(4)*ones(Ntime,1);  
            % IG =6; InitialGuess = H_range(5)*ones(Ntime,1);  
            % IG =7; indices = [3 5 1 2 4 2 3 5 1 3 2 4 5 2 1 4 5 3 1 5 4 2 3 1 5 2 3 1 5 2]; InitialGuess = H_range(indices);  % 1. Random variation of H bounded by the minimum and maximum values
            IG =8; [H_var_steps,] = H_variation_between_H_AND_0(Ntime); InitialGuess = 7600*H_var_steps;    % Rob's group H variation Peak amplitude was 7600 [A/m] with 50% duty cycle (60s On and 60s Off)


            pmin =  zeros(Ntime,1);
            pmax =  50000*ones(Ntime,1);
            tolx=1.e-9;
            tolfun=5.e-4;
            maxiter=400;

            % truthconstraint = infeasibility(stateconstraint,x0);
            %[popt,fval,exitflag,output] = solve(convprob,x0,'Options',myoptions, 'ConstraintDerivative', 'auto-reverse', 'ObjectiveDerivative', 'auto-reverse' )
            Fx = @(x) MIGHQuadMHT(x,NGauss,NumberUncertain,Nspecies,Ntime,GaussLegendre,ObjectiveType,deltat);
            [designopt,fval,exitflag,output,lambda,grad,hessian] ...
                =fmincon(Fx, InitialGuess ,[],[],[],[],pmin,pmax,[],myoptions);

            % handle = figure(5);

            %% optparams.FaList = reshape(designopt(:),size(params.FaList ));
            %% refparams.FaList = reshape(designopt(:),size(params.FaList ));
            %% popt.FaList      = reshape(designopt(:),size(params.FaList ));
    end
    % save convergence history
    set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman', 'FontSize', 12);
    savename1 = sprintf('IG%dNU%d_Fe3O4NP_ConvHistory',IG,NumberUncertain)
    saveas(gcf,savename1,'png')
    saveas(gcf,savename1,'fig')
    % saveas(handle,sprintf('historyNG%dNu%d%s%sSNR%02d%s',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule ),'png')
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
    %% plot initial and optimum H Vecotr

    figure('Position',[50,50,400,275])
    % subplot(1,2,1)
    steps = (1:1:Ntime)';
    stairs(InitialGuess,'r', 'LineWidth', 2);  % Using stairs for a more accurate representation
    hold on
    stairs(designopt, ':ob','LineWidth', 2);  % Using stairs for a more accurate representation
    hold off
    legend('Initial guess', 'Optimized',Location='best')
    % title('Initial and Optimum H variation');
    xlabel('Time Step Number');
    ylabel('Magnetic field amplitude, (A/m) ');
    ylim([0, 50000]);  % Set y-axis limits for better visibility
    xlim([0, Ntime]);   % Set x-axis limits, starting from 1
    grid on;
    set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman', 'FontSize', 12);
    savename2 = sprintf('IG%dNU%d_Fe3O4NP_OptH',IG,NumberUncertain)
    saveas(gcf,savename2,'png')
    saveas(gcf,savename2,'fig')
    
    figure('Position',[50,50,400,275])
    % subplot(1,2,2)
    stairs(grad,'r', 'LineWidth', 2);  % Using stairs for a more accurate representation
    xlabel('Time Step Number');
    ylabel('Gradient');
    % title('Gradient variation');
    cyl = ylim;
    ylim([cyl(1),0]);  % Set y-axis limits for better visibility
    xlim([0, Ntime]);   % Set x-axis limits, starting from 1
    set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman', 'FontSize', 12);
    savename3 = sprintf('IG%dNU%d_Fe3O4NP_Grad',IG,NumberUncertain)
    saveas(gcf,savename3,'png')
    saveas(gcf,savename3,'fig')
    % Htime =designopt;
    % %%% Time discretization: Implicit method is employed
    % dt = 1; % time-step used in for forward runs
    % t = 0:length(Htime)*deltat;     % Discrete times
    % TS = length(t);     % Number of time steps
    %
    % %%% Assigning the PieceWise Stepped H in deltat at respective times at dt resolution
    % for i = 1: length(Htime)
    %     startIndex = (i - 1) * (deltat)/dt + 1; % Calculate the starting index for the current iteration
    %     endIndex = min(startIndex + (deltat)/dt - 1, TS-1); % Calculate the ending index for the current iteration
    %     HVector(startIndex:endIndex) = Htime(i); % Amplitude vector in time at a resolution of dt assigned from deltat resolution
    % end
    % figure()
    % plot(t,[HVector(1),HVector])
    % xlim([0,t(end)])
    % xticks([0:t(end)/4:t(end)]);
    % xlabel('Time [s]')
    % ylabel('H Amplitude, [A/m]');
    % title('Optimized Temporal Magnetic Field');
end



%% evaluate MI
function MIobjfun =MIGHQuadMHT(hOpt,NGauss,NumberUncertain,Nspecies,Ntime,GaussLegendre,ObjectiveType,deltat)
%% Thermal conductivity
kmean = [ 0.492 ]; % [W/m/K]
kstdd = [ 0.076 ]; % [W/m/K]
klb   = [ 0.300 ]; % [W/m/K]
kub   = [ 0.570 ]; % [W/m/K]
%% perfusion
wmean = [ 0.00682 ]; % [1/s]
wstdd = [ 0.00387 ]; % [1/s]
wlb   = [ 0.00111 ]; % [1/s]
wub   = [ 0.01390 ]; % [1/s]
%% rho * specific heat
crhomean = [ 3589229 ];       % [J/m^3/K]
crhostdd = [ 0417740 ];       % [J/m^3/K]
crholb   = [ 2505204 ];       % [J/m^3/K]
crhoub   = [ 4056000 ];       % [J/m^3/K]

%% Anisotropy Constant

kmmean = [ 31500 ];     % [J/m^3]
kmstdd = [ 4725 ];      % [J/m^3]
kmlb   = [ 20000 ];      % [J/m^3]
kmub   = [ 53000 ];     % [J/m^3]

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
        %%% uncertain variable is K anisotropy constant
        if(GaussLegendre)
            [x,xn,xm,w,wn]=GaussLegendreNDGauss(NGauss,kmlb,kmub);
        else
            [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,kmmean,kmstdd);
        end
        kmqp   = xn{1}(:);

        % %%% uncertain variable is w blood perfusion rate
        % if(GaussLegendre)
        %     [x,xn,xm,w,wn]=GaussLegendreNDGauss(NGauss,wlb,wub);
        % else
        %     [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,wmean,wstdd);
        % end
        % wqp   = xn{1}(:);

        % evaluate function at each quadrature point
        lqp=length(xn{1}(:));
        sumstatevariable      = zeros(Nspecies,lqp);
        for iqp = 1:lqp
            sumstatevariable(1,iqp) =  pennesmht(kmean,wmean,crhomean,kmqp(iqp),hOpt,deltat);
            % sumstatevariable(1,iqp) =  pennesmht(kmean,wqp(iqp),crhomean,kmmean,hOpt,deltat);

        end
    case(4)
        if(GaussLegendre)
            [x,xn,xm,w,wn]=GaussLegendreNDGauss(NGauss,[klb; wlb; crholb; kmlb],[kub; wub; crhoub; kmub]);
        else
            [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,[kmean; wmean; crhomean; kmmean],[kstdd; wstdd; crhostdd; kmstdd]);
        end
        kqp    = xn{1}(:);
        wqp    = xn{2}(:);
        crhoqp = xn{3}(:);
        kmqp  = xn{4}(:);
        % evaluate function at each quadrature point
        lqp=length(xn{1}(:));
        sumstatevariable      = zeros(Nspecies,lqp);
        for iqp = 1:lqp
            sumstatevariable(1,iqp) =  pennesmht(kqp(iqp),wqp(iqp),crhoqp(iqp),kmqp(iqp),hOpt,deltat);
        end
end
%disp('build objective function')
expandvar  = ones(1,lqp);

switch (ObjectiveType)
    case('TotalSignal')
        diffsumm =sumstatevariable(1,:)' * expandvar   - expandvar' * sumstatevariable(1,:);
        negHz = 0;
        for jjj=1:lqp2
            znu=xn2{1}(jjj);
            % note that the sqrt(pi)^N+1  and 2^N factors from the integration over priors is included in the quadrature routines.
            % this makes is easier to switch between Gaussian and Uniform RV
            negHz = negHz + wn2(jjj) * (wn(:)' * log(exp(-(znu + diffsumm).^2/2/signu^2 - log(signu) -log(2*pi)/2) * wn(:)));
        end
end
% MI = H(z) - H(z|P)
%  H(z|P)  constant ==> max MI = max H(z) = min -H(z)
MIobjfun = negHz;
% MIobjfun = [kmqp',sumstatevariable,negHz];
end



%% Rob's Group H variation
function [rob_group_H,steps] = H_variation_between_H_AND_0(Ntime)
% Generate x values with a frequency of 1, starting from 1
steps = 1:1:Ntime;

% Calculate y values using the modified rectangle_wave function
rob_group_H = rectangle_wave(steps);
% Plot H variation
% figure;
% stairs(steps, rob_group_H, 'LineWidth', 2);  % Using stairs for a more accurate representation
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

%% Pennes Solver
function tempqoi = pennesmht(k_t,w_t,rhocp_t,Km,Htime,deltat)
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
t = 0:length(Htime)*deltat;     % Discrete times
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
T_SS = steady_state_temperature(T_b, htc, T_amb, Qm_t, rho_b, cp_b, w_t, k_t, r);
T_initial = T_SS;     % [°C] Initial Condition Temperature
% T_initial = 37;     % [°C] Initial Condition Temperature
alpha_t = k_t/(rhocp_t); % Tissue Thermal diffusivity [m^2/s]
f = 1.20e5;    % [Hz] Magnetic Field Frequency f = 163 kHz

% %%% This section has particle magnetic and physical properties for Fe nanoparticles being used in in-vivo experiments
% %%% Magnetic Nanoparticle Physical and Magnetic Parameters MNP Used: Fe3O4
% d_mnp = 60e-9;                 % [m] MNP Diameter d = 60 nm
% delta_l = 0.1*d_mnp;                 % [m] Liquid Layer Thickness on the Hard Solid MNP delta = 1 nm
% dh_mnp = d_mnp+2*delta_l;       % [m] MNP Hydrodynamic Diameter
% %%% domain magnetization
% % Mass specific magnetization of MNPs used in mice experiments is 150 emu/g
% % conversion of emu/g to [A/m]
% % 1 emu/g = 1Am2/kg
% % Am2/kg to A/m obtained by multiplying with density in kg/m3
% rho_mnp = 7874;        % [kg/m3] MNP Density
% Md_measured = 150;                  % [emu/g]
% Md = Md_measured*rho_mnp;         % [A/m] Domain Magnetization

%%% This section has particle magnetic and physical properties for Fe3O4 nanoparticles just to test
%%% Magnetic Nanoparticle Physical and Magnetic Parameters MNP Used: Fe3O4
d_mnp = 16e-9;                 % [m] MNP Diameter d = 60 nm
delta_l = 0.1*d_mnp;                 % [m] Liquid Layer Thickness on the Hard Solid MNP delta = 1 nm
dh_mnp = d_mnp+2*delta_l;       % [m] MNP Hydrodynamic Diameter
%%% domain magnetization
% Mass specific magnetization of MNPs used in mice experiments is 150 emu/g
% conversion of emu/g to [A/m]
% 1 emu/g = 1Am2/kg
% Am2/kg to A/m obtained by multiplying with density in kg/m3
rho_mnp = 5180;        % [kg/m3] MNP Density
Md = 446*1000;         % [A/m] Domain Magnetization


%%% MNP Dose and Magnetic Fluid Parameters
mnp_tumor_conc = 4; % %[kg/m^3] kg of MNP/m3 of tumor = 1mg/cm3 [ug/ml] micrograms of MNP/milliliter of water
mf_vol = 1; %[ml] microliters of nanofluid injected
mf_vol2 = mf_vol*1e-6; % [m3] of nanofluid
mnp_mass_kg = mnp_tumor_conc*vol_t;
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
gamma = Km*mnp_svol/(kB*T_bK);                    % [-] An intermediate parameters used in further calculations
tauB = (3*mu_cf*mnp_hvol)/(kB*T_bK);                % [s] Brownian Relaxation Time of MNP
tauN = tau0*(sqrt(pi)/2)*(1/sqrt(gamma))*exp(gamma) ;   % [s] Neel Relaxation Time of MNP
% if isinf(tauN)
%     error(['tauN is infinity and gamma = ', num2str(gamma)])
% end
One_by_tauE = 1/tauB+1/tauN;                       % [s] Effective Relaxation Time of MNP
tauE = 1/One_by_tauE;

%%% Implicit Finite Difference Scheme using the Cauchy Boundary Condition on the outer boundary
T = zeros(N, TS);
T(:,1) = T_initial; % Initial Condition
% Coefficient Matrix Three-Diagonal Generation
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
% Solver Time Iterations
for n = 2:TS
    H = HVector(n-1); %
    if H == 0
        SAR = 0;
    else
        zeta = mu0*Md*H*mnp_svol/(kB*T_bK);                 % [-] Langevin parameter
        X_i = mu0*(Md^2)*mnp_vol_frac*mnp_svol/(3*kB*T_bK); % [-] Initial Magnetic Susceptibility
        X_0 = 3*X_i*(coth(zeta)-1/zeta)/zeta;               % [-] Equilibrium Magnetic Susceptibility
        X_L = X_0*omega*tauE/(1+(omega*tauE)^2);            % [-] Loss Component of Magnetic Susceptibility
        P = pi*mu0*f*X_L*H^2;                               % [W/m3] Heat Generation rate of MNPs in MF
        SAR = P/(rho_mnp*mnp_vol_frac);                     % [W/kg] MNP Specific Absorption Rate
        SAR_grams = SAR/1000;                              % [W/g] MNP Specific Absorption Rate
    end
    Q_MNP = alpha_CF*mnp_tumor_conc*SAR;  % [W/m3] Heat Generation by MNP in Tissue

    %%% Assign Q_MNP to the center region only
    q_mnp = Q_MNP*t_loc;        % MNP heat generation [W/m^3] within the tumor only
    q_mnp_time(:,n-1) = q_mnp;
    Force = zeros(N,1); % Force vector
    Force(1:N-1)  = T(1:N-1,n-1) + (rho_b*cp_b.*w_t*dt*T_b + Qm_t*dt+ q_mnp(1:N-1)*dt)/(rhocp_t);
    Force(N) = htc*T_amb;
    T(:,n) = thomas_algorithm(Lower, Main, Upper, Force); % Spatio-Temporal Temperature
end


% Spatial integration (integrating over radial locations for each time point)
spatial_integral = trapz(r, T, 1);  % Integrate over the radial direction

% Time integration (integrating the result over time)
total_integral = trapz(t, spatial_integral);  % Integrate over time

% T_sum_time = sum(T,2);      % Summation of Temperature over all Times at each r location
% G = trapz(r, T_sum_time);     % Integration of Temperature over the domain in 3D
tempqoi = total_integral;
%%% Plot the results
% q_mnp_t = [q_mnp_time(:,1),q_mnp_time];
% [RR,TT] = meshgrid(r,t);
% figure('Position', [40, 40, 800, 600])
% sgtitle('Fe_3O_4 NP Bioheat Solution')
% subplot(2,2,1)
% plot((1:length(HVector))*dt,HVector)
% xlim([0,t(end)])
% xticks([0:t(end)/4:t(end)]);
% xlabel('Time [s]')
% ylabel('H Amplitude, [A/m]');
% title('Magnetic Field');
% 
% subplot(2,2,2)
% surf(RR,TT,q_mnp_t')
% xlabel('Radial Distance [m]')
% ylabel('Time [s]');
% zlabel('Qmnp [W/m^3]');
% title('Spatio-Temporal Heat Source');
% xlim([0,R])
% xticks([0,RT/2,RT,R]);
% ylim([0,t(end)])
% yticks([0:t(end)/4:t(end)]);
% grid on;
% 
% subplot(2,2,3)
% figure(505)
% plot_radius = [0,RT,R];
% legend_String = string(plot_radius)+[" m Tumor Center"," m Tumor Edge"," m Outer Boundary"];
% plot_rad_idx = (plot_radius/dr)+1;
% plot(t,T(plot_rad_idx,:), 'LineWidth',2)
% xlabel('Time, t [s]');
% xlim([0,t(end)])
% xticks([0:t(end)/4:t(end)]);
% cylims = ylim;
% ylim ([floor(cylims(1)),ceil(cylims(2))])
% ylabel('Temperature, T [°C]');
% legend(legend_String, Location='best')
% title('Temperature Elevations');
% grid on;
% 
% subplot(2,2,4)
% plot_time = [0,60,120,t(end)];
% legend_String = string(plot_time)+[" s"," s"," s"," s"];
% plot_ind = (plot_time/dt)+1;
% plot(r,T(:,plot_ind),'LineWidth',2)
% xlabel('Radial distance, r [m]');
% xlim([0,R])
% xticks(0:R/5:R);
% ylabel('Temperature, T [°C]');
% cylims = ylim;
% ylim ([floor(cylims(1)),ceil(cylims(2))])
% legend(legend_String, Location='best')
% title('Spatial Temperature profile');
% grid on;
% alltextObjects = findall(gcf,'-property','FontName');
% set(alltextObjects, 'FontName','Times new Roman', 'FontSize',12)
% saveas(gca,['Fe3O4SolutionForIG_',num2str(IG)],'png')
%%% Tri-Diagonal Matrix Alogorithm for Linear System of Equation Solver
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

%%% Steady-State Temperature Profiles without the external heating, to calculate initial condictions for the transient heat transfer model
    function T_SS = steady_state_temperature(T_b, htc, T_amb, Qm_t, rho_b, cp_b, w_t, k_t, r)
        % Inputs:
        % T_b: Blood temperature (°C)
        % htc: Heat transfer coefficient (W/m^2K)
        % T_amb: Ambient temperature (°C)
        % Qm_t: Metabolic heat generation rate (W/m^3)
        % rho_b: Blood density (kg/m^3)
        % cp_b: Blood specific heat capacity (J/kgK)
        % w_t: Blood perfusion rate (1/s)
        % k_t: Tissue thermal conductivity (W/mK)
        % R: Radius of the domain (m)
        % N: Number of radial grid points

        % Discretization parameters
        delta_r = r(2)-r(1);        % Radial step size
        Num_ele = length(r); % Radial positions

        % Initialize system matrix A and right-hand side vector b
        A = zeros(Num_ele, Num_ele);
        b = zeros(Num_ele, 1);

        % Coefficients for blood perfusion and metabolic heat generation
        perfusion_term = w_t * rho_b * cp_b;

        % Fill the matrix A and vector b for internal nodes (i = 2 to N-1)
        for idx = 2:Num_ele-1
            r_i = r(idx);

            % Coefficients for the discretized equation
            coef_T_ip1 = (k_t / delta_r^2) + (k_t / (r_i *  delta_r)); % T_{i+1}
            coef_T_im1 = (k_t / delta_r^2) - (k_t / (r_i *  delta_r)); % T_{i-1}
            coef_T_i = -(2 * k_t / delta_r^2 + perfusion_term);       % T_i

            % Populate the matrix A and vector b
            A(idx, idx+1) = coef_T_ip1;      % T_{i+1}
            A(idx, idx-1) = coef_T_im1;      % T_{i-1}
            A(idx, idx) = coef_T_i;          % T_i
            b(idx) = -perfusion_term * T_b - Qm_t;
        end

        % Boundary condition at r = 0 (Symmetry: T_1 = T_0)
        A(1, 1) = 1;
        A(1, 2) = -1;
        b(1) = 0;

        % Boundary condition at r = R (Convective boundary condition)
        A(Num_ele, Num_ele-1) = k_t / delta_r;
        A(Num_ele, Num_ele) = -(k_t / delta_r + htc);
        b(Num_ele) = -htc * T_amb;

        % Solve the system of linear equations A*T = b
        T_SS = A \ b;

        % Plot the temperature profile
        % figure;
        % plot(r, T_SS, '-o');
        % xlabel('Radial Position (m)');
        % ylabel('Temperature (°C)');
        % title('Temperature Profile in the Tissue');
        % grid on;
    end
end