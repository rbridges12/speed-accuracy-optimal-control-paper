function result = optimization_6muscles(N, wM_std, pos_conf_95, vel_conf_95, k_u, k_t, X_init, u_init, hot_start, dX_init, Pmat_init, EE_target)

    import casadi.*

    % Set-up structure with data that specifies model
    auxdata = initializeModelParameters();

    % Additional simulation settings
    % T = 0.8;
    auxdata.N = N; % number of discretized nodes
    nStates = 4; auxdata.nStates = nStates; % [q1, q2, qdot1, qdot2]
    nControls = auxdata.nMuscles;

    %%%% Define CasADi functions - for reasons of efficiency and to compute sensitivities (jacobians) of functions
    functions = generateFunctions(auxdata);

    opti = casadi.Opti(); % Create opti instance

    % create optimization variables and provide initial guesses
    T_guess = 0.8;
    X_guess = zeros(nStates,N+1);
    X = opti.variable(nStates,N+1);
    opti.set_initial(X, X_guess);
    u = opti.variable(nControls, N+1);
    opti.set_initial(u, 0.01);
    % if (hot_start)
    %     opti.set_initial(u, u_init);
    % end
    M = opti.variable(nStates,nStates*N);
    opti.set_initial(M, 0.01);
    T = opti.variable(); % Duration
    opti.set_initial(T, T_guess);

    dt = T/N;
    % time = 0:dt:T;
    time = linspace(0,T,N+1);
    wM = (wM_std*ones(auxdata.nMotorNoises,1)).^2/dt; % Motor noise: go from std of continuous noise source to variance of discrete sample
    sigma_w = diag(wM);

    Pmat_i = Pmat_init;
    for i = 1:N
        X_i = X(:,i);
        X_i_plus = X(:,i+1);
        u_i = u(:,i);
        u_i_plus = u(:,i+1);
        
        % dynamics constraints using trapezoidal integration
        dX_i = functions.f_forwardMusculoskeletalDynamics(X_i, u_i, 0);
        dX_i_plus = functions.f_forwardMusculoskeletalDynamics(X_i_plus, u_i_plus, 0);
        opti.subject_to(functions.f_G_Trapezoidal(X_i,X_i_plus,dX_i,dX_i_plus, dt)*1e3 == 0);
        
        M_i = M(:,(i-1)*nStates + 1:i*nStates);
        
        DdX_DX_i = functions.f_DdX_DX(X_i,u_i,wM);
        DdZ_DX_i = functions.f_DdX_DX(X_i_plus,u_i_plus,wM);
        DdX_Dw_i = functions.f_DdX_Dw(X_i,u_i,wM);
        DdZ_Dw_i = functions.f_DdX_Dw(X_i_plus,u_i_plus,wM);
        % DdX_DX_i = functions.f_DdX_DX(X_i,u_i,0);
        % DdZ_DX_i = functions.f_DdX_DX(X_i_plus,u_i_plus,0);
        % DdX_Dw_i = functions.f_DdX_Dw(X_i,u_i,0);
        % DdZ_Dw_i = functions.f_DdX_Dw(X_i_plus,u_i_plus,0);
        
        DG_DX_i = functions.f_DG_DX(DdX_DX_i, dt);
        DG_DZ_i = functions.f_DG_DZ(DdZ_DX_i, dt);
        DG_DW_i = functions.f_DG_DW(DdX_Dw_i, DdZ_Dw_i, dt);
        
        opti.subject_to(M_i*DG_DZ_i - eye(nStates) == 0);
        
        Pmat_i = M_i*(DG_DX_i*Pmat_i*DG_DX_i' + DG_DW_i*sigma_w*DG_DW_i')*M_i';
    end

    % constrain time to be positive
    opti.subject_to(T > 0);

    % constrain initial and final conditions
    opti.subject_to(X(:,1) - X_init == 0);
    % opti.subject_to(functions.f_EEPos(X(1:2, 1)) - q_init == 0);
    % opti.subject_to(functions.f_EEPos(X(1:2, 1)) - EE_init == 0);
    dX_init_opt = functions.f_forwardMusculoskeletalDynamics(X(:, 1), u(:, 1), 0);
    if (hot_start)
        opti.subject_to(u(:,1) - u_init == 0);
    end
    opti.subject_to(dX_init_opt == dX_init); % Initial acceleration equals zero (activations need to fullfill requirement)

    % Reaching motion must end in the final reach position
    opti.subject_to(functions.f_EEPos(X(1:2,end)) - EE_target == 0);

    % Final joint velocity and acceleration equals zero (activations balanced)
    dX_end = functions.f_forwardMusculoskeletalDynamics(X(:,end),u(:,end),0);
    opti.subject_to(dX_end == 0);

    
    % Constrain end point position and velocity variance in x and y
    P_q_final = Pmat_i(1:2, 1:2);
    P_EEPos_final = functions.f_P_EEPos(X(1:2,end),P_q_final);
    P_q_qdot_final = Pmat_i;
    P_EEVel_final = functions.f_P_EEVel(X(1:2,end),X(3:4,end),P_q_qdot_final);
    pos_variance = (pos_conf_95/2)^2;
    vel_variance = (vel_conf_95/2)^2;
    opti.subject_to(P_EEPos_final(1,1) < pos_variance);
    opti.subject_to(P_EEPos_final(2,2) < pos_variance);
    opti.subject_to(P_EEVel_final(1,1) < vel_variance);
    opti.subject_to(P_EEVel_final(2,2) < vel_variance);

    % constrain activations and joint angle limits
    opti.subject_to(0.001 < u(:) < 1);
    opti.subject_to(0 < X(1,:) < 180);
    opti.subject_to(0 < X(2,:) < 180);

    % minimize weighted combination of energy usage (integral of activations) and duration
    opti.minimize(k_u*1e3*(sumsqr(u)/2)*dt + k_t*T);

    % Set solver options
    % optionssol.ipopt.nlp_scaling_method = 'gradient-based';
    optionssol.ipopt.linear_solver = 'mumps';
    optionssol.ipopt.tol = 1e-3;
    optionssol.ipopt.dual_inf_tol = 3e-4;
    optionssol.ipopt.constr_viol_tol = 1e-7;
    optionssol.ipopt.max_iter = 1000;
    optionssol.ipopt.hessian_approximation = 'limited-memory';
    opti.solver('ipopt',optionssol);

    sol = opti.solve();
    u_sol = sol.value(u);
    X_sol = sol.value(X);
    X_init_sol = X_sol(:,1);
    T_sol = sol.value(T);
    dt_sol = T_sol/N;
    wM_sol = (wM_std*ones(auxdata.nMotorNoises,1)).^2/dt_sol;
    sigma_w_sol = diag(wM_sol);
    auxdata.wM = wM_sol;
    auxdata.sigma_w = sigma_w_sol;

    % simulate forward dynamics with initial P to compute covariance matrices
    [X_sim, ~, EE_ref_sol, Pmat_sol] = forwardSim_no_fb(X_init_sol,Pmat_init,u_sol, T_sol,auxdata,functions);

    final_cost = sol.value(opti.f);

    % Compute end-effector position and velocity covariances
    for i = 1:N+1
        Pmat_sol_i = Pmat_sol(:,:,i);
        P_q_sol_i = Pmat_sol_i(1:2, 1:2);
        P_EEPos_mat_i = functions.f_P_EEPos(X_sim(1:2,i),P_q_sol_i);
        P_EEPos_sol(:,i) = full([P_EEPos_mat_i(1,1); P_EEPos_mat_i(2,1); P_EEPos_mat_i(2,2)]);
        P_qdot_sol_i = Pmat_sol_i;
        P_EEVel_mat_i = functions.f_P_EEVel(X_sim(1:2,i),X_sim(3:4,i),P_qdot_sol_i);
        P_EEVel_sol(:,i) = full([P_EEVel_mat_i(1,1); P_EEVel_mat_i(2,1); P_EEVel_mat_i(2,2)]);
    end

    EEPos_sol = EndEffectorPos(X_sim(1:2,:),auxdata)';
    EEVel_sol = EndEffectorVel(X_sim(1:2,:),X_sim(3:4,:),auxdata)';

    final_pos_cov = P_EEPos_sol(:,end);
    final_vel_cov = P_EEVel_sol(:,end);

    clear result;
    result.e_ff = u_sol';
    result.X = X_sim;
    result.q = X_sol(1:2,:)';
    result.qdot = X_sol(3:4,:)';
    result.M = sol.value(M);
    result.Pmat = Pmat_sol;
    result.time = 0:dt_sol:T_sol;
    result.auxdata = auxdata;
    result.functions = functions;
    result.EEPos = EEPos_sol;
    result.EEVel = EEVel_sol;
    result.P_EEPos = P_EEPos_sol;
    result.P_EEVel = P_EEVel_sol;
    result.EE_ref = EE_ref_sol;
    result.final_cost = final_cost;
    result.final_pos_cov = final_pos_cov;
    result.final_vel_cov = final_vel_cov;
    result.EE_target = EE_target;
    result.target_width = pos_conf_95;
    result.Sigma_w = sigma_w_sol;
    result.debug = opti.debug;
    result.stats = opti.stats;
end



