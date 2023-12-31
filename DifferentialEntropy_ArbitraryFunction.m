function [diff_ent] = DifferentialEntropy_ArbitraryFunction(mu_theta,cov_theta,dim_theta,dim_z,f_theta,cov_z,N_samples_theta,N_samples_z,k,N_trials)
% The function explores the calculation of an integral of the form 
%               H(z) = -\int_Z p(z) \log(p(z)) dz
% under the condition that p(z) is not itself well-defined, but instead 
% expressed as 
%               p(z) = \int_{\Theta} p(z|\theta) p(\theta) d\theta
% In this function, we'll assume that 
%   (1) p(theta) is Gaussian with mean \mu_theta and cov \Sigma_theta
%   (2) p(z|theta) is conditionally Gaussian with mean \mu_Z = f(\theta) for some
%           user defined function f, and cov \Sigma_z
% The shortcuts we use here are that:
%   (1) since theta is a standard normal, we generate those samples using randn. 
%   (2) since z | theta is normal, we can generate it using randn as well 
% Function Inputs:
%   - mu_theta: the mean of theta
%   - cov_theta: the cov matrix of theta
%   - dim_theta: the dimension of theta, a real Gaussian-dist. vector
%   - dim_z: the dimension of z, a real Gaussian-dist. vector
%   - @f_theta: a function for calculating the mean of z from theta
%   - cov_z: the cov matrix of z
%   - N_samples_theta: the number of samples for the inner loop of the MC
%   `   computation, which makes realizations theta_j
%   - N_samples_outer: the number of samples for the outer loop, which
%   makes realizations of z_i|theta_j

% Function Outputs:
%   - diff_ent: the estimated entropy of z 
if size(mu_theta,1) ~= dim_theta 
    error('Please specify mean vector of theta as dim_theta x 1 vector');
end
if ~all(size(cov_theta) == [dim_theta,dim_theta])
    error('Please specify cov mtx as (dim_theta x dim_theta) mtx'); 
end
if ~all(size(cov_z) == [dim_z, dim_z])
    error('Please specify cov matrix of z as (dim_z x dim_z)'); 
end


mu_theta = transpose(mu_theta); % resolves dimension with 'start_point'
burn_in_N = 1000;

diff_ent_measured = zeros(N_trials,1); 
for ii = 1:N_trials 
    
    kt = k'
    %% sample from p(theta)
    % theta ~ N(mu_theta, cov_theta)
    theta_i = cov_theta*(randn([dim_theta,N_samples_theta])) + mu_theta*ones([dim_theta, N_samples_theta]); 
    f_theta_i = f_theta(theta_i,k);

    %% sample from p(z|theta)
    theta_for_z = cov_theta*(randn([dim_theta,N_samples_z])) + mu_theta*ones([dim_theta, N_samples_z]); 
    theta_z_generation =   f_theta(theta_for_z,k);  
    z_given_theta = (cov_z*randn([dim_z,N_samples_z]) + theta_z_generation)';
    det_cov_z = det(cov_z); % pre-compute determinant and inverse 
    inv_cov_z = inv(cov_z); 

    %% Compute differential entropy 
    % H(z) = -\int_z p(z) log(p(z))
    %      
    % p(z) = E_{\Theta} [ p(z|\theta_i)] 
        
    % define computation of p(z)
    pdf_z = @(z) mean(compute_p_z(z, f_theta_i, det_cov_z, inv_cov_z)); % computes p(z) = E_{theta}[p(z | theta)]    
    
    % calculate E_{\Theta}[p(z|\theta_i)] 
    p_z = zeros(N_samples_z,1);
    for i = 1:size(z_given_theta,1)
        z_i = z_given_theta(i,:);
        p_z(i) = pdf_z(z_i);
    end
    log_pdf_values = log(p_z);
    diff_ent = -1*mean(log_pdf_values); % = E[-log(p(x))]
    diff_ent
    diff_ent_measured(ii) = diff_ent; 
%     toc
end
end

%% Given realizations of theta, compute p(z|theta)
% inputs:
%   - z_0: single realization of z, used to compute p(z|theta_i)
%   - theta_realizations: realizations of theta, i.e. theta_i
%   - f_mu_z_theta: the function for computing \mu_z from theta
%   - sigma_z: the covar of z
% outputs
%   - p(z) for given z
% description
%   - goal: p(z) = \int_{\Theta} p(z|theta) p(theta) d\theta
%                \approx \frac{1}{N} \sum_i p(z | theta_i)
%   1. for each theta_i, compute \mu_z = f(theta_i) = f_mu_z_theta(theta_i)
%   2. this sets a distribution z|theta_i ~ N(\mu_z, sigma_z)
%   3. use mvnpdf with those params to eval p(z|theta_i)
%   4. take the average of this over many theta_i
%   5. this average corresponds to an estimate for p(z) 
function p_z = compute_p_z(z_0, f_theta, det_sigma_z, inv_sigma_z)
    N = size(f_theta,2);
    dim = size(z_0,2); 
    p_z_given_theta = zeros(N,1); 
    f_theta = f_theta.'; 
    for i = 1:N
        p_z_given_theta(i) = (2*pi)^(-dim/2)*det_sigma_z^(-0.5)*exp(-0.5*(z_0 - f_theta(i,:))*inv_sigma_z*(z_0 - f_theta(i,:)).');
        
        % deprecated approach which recomputed inverse and determinant each time 
        %mvnpdf(z_0,f_theta(i,:),sigma_z); % the objective function is p(z | theta) ~ \mathcal{N}(f(theta), Sigma_z) which is to say that
                                                              % the mean of z is conditional on theta through f(theta)
    end
    p_z = mean(p_z_given_theta); % p(z) \approx \frac{1}{N} \sum_i p(z | \theta_i)
end

