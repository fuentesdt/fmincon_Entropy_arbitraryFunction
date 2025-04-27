function temperature_profile = predict_temperature_profile(T_b, htc, T_amb, Qm_t, rho_b, cp_b, w_t, k_t, R, r)
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
        % N: Number of discretized radial points

        % Discretization
        del_r = r(2)-r(1);       % Radial step size
        Num_el = length(r);
        % Preallocate temperature array
        T_Steady = zeros(Num_el, 1);

        % Coefficients for the Pennes bioheat equation
        A = @(idx) k_t * (r(idx+1)^2) / (r(idx)^2 * del_r^2);  % Coefficient for T_{i+1}
        B = @(idx) k_t * (r(idx-1)^2) / (r(idx)^2 * del_r^2);  % Coefficient for T_{i-1}
        C = @(idx) w_t * rho_b * cp_b;                  % Blood perfusion term

        % System of equations matrix
        A_mat = zeros(Num_el, Num_el);
        b_vec = zeros(Num_el, 1);

        % Boundary condition at r = 0 (Symmetry condition: T_1 = T_0)
        A_mat(1, 1) = 1;
        A_mat(1, 2) = -1;

        % Boundary condition at r = R (Convective boundary condition)
        A_mat(Num_el, Num_el-1) = k_t / del_r;
        A_mat(Num_el, Num_el) = -(k_t / del_r + htc);
        b_vec(Num_el) = -htc * T_amb;

        % Fill the interior points (for i = 2 to N-1)
        for idx = 2:Num_el-1
            A_mat(idx, idx+1) = A(idx);
            A_mat(idx, idx-1) = B(idx);
            A_mat(idx, idx) = -(A(idx) + B(idx) + C(idx));
            b_vec(idx) = -C(idx) * T_b - Qm_t;
        end

        % Solve the linear system A_mat * T = b_vec
        T_Steady = A_mat \ b_vec;

        % Return the temperature profile
        temperature_profile = T_Steady;

        % Plot the temperature profile
        figure;
        plot(r, T_Steady, '-o');
        xlabel('Radial Position (m)');
        ylabel('Temperature (°C)');
        title('Temperature Profile in the Tissue');
        grid on;
    end