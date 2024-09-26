function q_sol = ik_opt(EE_target)

    import casadi.*

    % Set-up structure with data that specifies model
    auxdata = initializeModelParameters();

    % Additional simulation settings
    % T = 0.8;
    auxdata.N = 80; % number of discretized nodes
    nStates = 4; auxdata.nStates = nStates; % [q1, q2, qdot1, qdot2]
    nControls = auxdata.nMuscles;

    %%%% Define CasADi functions - for reasons of efficiency and to compute sensitivities (jacobians) of functions
    functions = generateFunctions(auxdata);

    opti = casadi.Opti(); % Create opti instance

    % create optimization variables and provide initial guesses
    q = opti.variable(2,1);
    opti.set_initial(q, 0);

    % opti.subject_to(functions.f_EEPos(q) - EE_target == 0);
    % opti.subject_to(0 < q(1) < 180);
    % opti.subject_to(0 < q(2) < 180);
    opti.subject_to(0 < q(1) < pi);
    opti.subject_to(0 < q(2) < pi);
    % opti.subject_to(q(2) < pi);
    % opti.subject_to(q(2) == 1);

    % minimize weighted combination of energy usage (integral of activations) and duration
    opti.minimize(100*norm(functions.f_EEPos(q) - EE_target));

    % Set solver options
    % optionssol.ipopt.nlp_scaling_method = 'gradient-based';
    % optionssol.ipopt.print_level = 0;
    optionssol.ipopt.linear_solver = 'mumps';
    optionssol.ipopt.tol = 1e-3;
    optionssol.ipopt.dual_inf_tol = 3e-4;
    optionssol.ipopt.constr_viol_tol = 1e-7;
    optionssol.ipopt.max_iter = 3000;
    optionssol.ipopt.hessian_approximation = 'limited-memory';
    opti.solver('ipopt',optionssol);

    sol = opti.solve();
    q_sol = sol.value(q);
end



