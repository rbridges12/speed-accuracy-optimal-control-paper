function result = optimization_6muscles(forceField,wM_std)

    import casadi.*

    % Set-up structure with data that specifies model
    auxdata = initializeModelParameters();

    % Additional simulation settings
    T = auxdata.T;
    dt = 0.01; auxdata.dt = dt; % time step
    N = round(T/dt); auxdata.N = N; % number of discretized nodes
    time = 0:dt:T; auxdata.time = time;
    nStates = 4; auxdata.nStates = nStates; % [q1, q2, qdot1, qdot2]
    wM = (wM_std*ones(2,1)).^2/dt; auxdata.wM = wM; % Motor noise: go from std of continuous noise source to variance of discrete sample

    sigma_w = diag(wM);
    auxdata.sigma_w = sigma_w;

    auxdata.forceField = forceField;

    %%%% Define CasADi functions - for reasons of efficiency and to compute sensitivities (jacobians) of functions
    functions = generateFunctions(auxdata);

    opti = casadi.Opti(); % Create opti instance

    % given shoulder angles, solve for initial and final elbow angles such
    % that the EE x position is 1 at the start and finish
    fsolve_options = optimoptions('fsolve','FiniteDifferenceType','central','StepTolerance',1e-10,'OptimalityTolerance',1e-10);
    shoulder_pos_init = 20*pi/180;
    shoulder_pos_final = 55*pi/180;
    f = @(x)get_px(x,auxdata,shoulder_pos_init);
    elbow_pos_init = fsolve(f,ones,fsolve_options);
    initial_pos = [shoulder_pos_init; elbow_pos_init];

    f = @(x)get_px(x,auxdata,shoulder_pos_final);
    elbow_pos_final = fsolve(f,ones,fsolve_options);
    final_pos = [shoulder_pos_final; elbow_pos_final];
    EE_final = EndEffectorPos(final_pos,auxdata);

    % create optimization variables and provide initial guesses
    X_guess = zeros(nStates,N+1);
    X_guess(1,:) = interp1([0 T], [initial_pos(1) final_pos(1)],time);
    X_guess(2,:) = interp1([0 T], [initial_pos(2) final_pos(2)],time);
    X = opti.variable(nStates,N+1);
    opti.set_initial(X, X_guess);
    u = opti.variable(6,N+1);
    opti.set_initial(u, 0.01);
    M = opti.variable(nStates,nStates*N);
    opti.set_initial(M, 0.01);
    Pmat_init = diag([1e-4; 1e-4; 1e-7; 1e-7]);

    Pmat_i = Pmat_init;
    for i = 1:N
        X_i = X(:,i);
        X_i_plus = X(:,i+1);
        u_i = u(:,i);
        u_i_plus = u(:,i+1);
        
        % dynamics constraints using trapezoidal integration
        dX_i = functions.f_forwardMusculoskeletalDynamics(X_i, u_i, 0);
        dX_i_plus = functions.f_forwardMusculoskeletalDynamics(X_i_plus, u_i_plus, 0);
        opti.subject_to(functions.f_G_Trapezoidal(X_i,X_i_plus,dX_i,dX_i_plus)*1e3 == 0);
        
        M_i = M(:,(i-1)*nStates + 1:i*nStates);
        
        DdX_DX_i = functions.f_DdX_DX(X_i,u_i,wM);
        DdZ_DX_i = functions.f_DdX_DX(X_i_plus,u_i_plus,wM);
        DdX_Dw_i = functions.f_DdX_Dw(X_i,u_i,wM);
        DdZ_Dw_i = functions.f_DdX_Dw(X_i_plus,u_i_plus,wM);
        
        DG_DX_i = functions.f_DG_DX(DdX_DX_i);
        DG_DZ_i = functions.f_DG_DZ(DdZ_DX_i);
        DG_DW_i = functions.f_DG_DW(DdX_Dw_i, DdZ_Dw_i);
        
        opti.subject_to(M_i*DG_DZ_i - eye(nStates) == 0);
        
        Pmat_i = M_i*(DG_DX_i*Pmat_i*DG_DX_i' + DG_DW_i*sigma_w*DG_DW_i')*M_i';
    end

    % constrain initial and final conditions
    opti.subject_to(X(:,1) - [initial_pos; 0; 0] == 0);
    dX_init = functions.f_forwardMusculoskeletalDynamics(X(:,1),u(:,1),0);
    opti.subject_to(dX_init([3:4]) == 0); % Initial acceleration equals zero (activations need to fullfill requirement)

    % Reaching motion must end in the final reach position with zero angular joint velocity
    opti.subject_to(functions.f_EEPos(X(1:2,end)) - EE_final == 0);
    opti.subject_to(X(3:4,end) == [0; 0]);

    % Final acceleration equals zero (activations balanced)
    dX_end = functions.f_forwardMusculoskeletalDynamics(X(:,end),u(:,end),0);
    opti.subject_to(dX_end(3:4) == 0);

    
    % Constrain end point position and velocity variance in x and y
    P_q_final = Pmat_i(1:2, 1:2);
    P_EEPos_final = functions.f_P_EEPos(X(1:2,end),P_q_final);
    P_q_qdot_final = Pmat_i;
    P_EEVel_final = functions.f_P_EEVel(X(1:2,end),X(3:4,end),P_q_qdot_final);
    pos_dev = 0.04^2;
    vel_dev = 0.5^2;
    opti.subject_to(P_EEPos_final(1,1) < pos_dev);
    opti.subject_to(P_EEPos_final(2,2) < pos_dev);
    opti.subject_to(P_EEVel_final(1,1) < vel_dev);
    opti.subject_to(P_EEVel_final(2,2) < vel_dev);

    % constrain activations and joint angle limits
    opti.subject_to(0.001 < u(:) < 1);
    opti.subject_to(0 < X(1,:) < 180);
    opti.subject_to(0 < X(2,:) < 180);
    % opti.subject_to((0.001 < u(:)) && (u(:) < 1));
    % opti.subject_to((0 < X(1,:)) && (X(1,:) < 180));
    % opti.subject_to((0 < X(2,:)) && (X(2,:) < 180));

    % minimize sum of squared activations
    opti.minimize(1e3*(sumsqr(u)/2)*dt);

    % Set solver options
    % optionssol.ipopt.nlp_scaling_method = 'gradient-based';
    optionssol.ipopt.linear_solver = 'mumps';
    optionssol.ipopt.tol = 1e-3;
    optionssol.ipopt.dual_inf_tol = 3e-4;
    optionssol.ipopt.constr_viol_tol = 1e-7;
    optionssol.ipopt.max_iter = 10000;
    optionssol.ipopt.hessian_approximation = 'limited-memory';
    opti.solver('ipopt',optionssol);

    sol = opti.solve();
    u_sol = sol.value(u);
    X_sol = sol.value(X);
    X_init_sol = X_sol(:,1);
    [X_sim, ~, EE_ref_sol, Pmat_sol] = forwardSim_no_fb(X_init_sol,Pmat_init,u_sol,auxdata,functions);

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

    clear result;
    result.e_ff = u_sol';
    result.X = X_sim;
    result.q = X_sol(1:2,:)';
    result.qdot = X_sol(3:4,:)';
    result.M = sol.value(M);
    result.Pmat = Pmat_sol;
    result.time = 0:dt:T;
    result.auxdata = auxdata;
    result.EEPos = EEPos_sol;
    result.EEVel = EEVel_sol;
    result.P_EEPos = P_EEPos_sol;
    result.P_EEVel = P_EEVel_sol;
    result.EE_ref = EE_ref_sol;
end



