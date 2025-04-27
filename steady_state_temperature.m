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
