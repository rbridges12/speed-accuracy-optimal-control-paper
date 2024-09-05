function result = nonlinear_mpc(N, motor_noise_stddev, target_radius, target_vel_accuracy, k_u, k_t, X_init, target_pos)
    Tsim = 0.1;
    max_iter_cold_start = 2000;
    max_rollouts = 50;
    P_0 = diag([1e-4; 1e-4; 1e-7; 1e-7]);

    P_init = P_0;
    current_time = 0;
    dX_init = zeros(4, 1);
    u_init = zeros(6, 1);
    hot_start = false;
    hs_X = [];
    hs_u = [];
    hs_M = [];
    hs_T = 0;
    max_iter = max_iter_cold_start;
    ts = [];
    ts_plans = [];
    X_plans = [];
    lengths = [];
    x_traj = [];
    u_traj = [];
    P_traj = [];
    ee_traj = [];
    P_EEPos = [];
    P_EEVel = [];

    for i = 1:max_rollouts
        result = optimization_6muscles(N, motor_noise_stddev, target_radius, target_vel_accuracy, k_u, k_t, X_init, u_init, dX_init, P_init, target_pos, hot_start, hs_X, hs_u, hs_M, hs_T, max_iter);

        dt = result.time(2) - result.time(1);
        Nsim = min(ceil(Tsim / dt), N);
        u_sim = result.e_ff';
        u_sim = u_sim(:, 1:Nsim+1);
        [X_sim, ~, EE_ref_sol, Pmat_sol] = forwardSim_ode(result.X(:,1), result.Pmat(:,:,1) ,result.e_ff', Nsim*dt, Nsim, motor_noise_stddev, result.auxdata, result.functions);

        ts_planned = result.time + current_time;
        ts_plans = [ts_plans; ts_planned];
        X_plans = cat(3, X_plans, result.X);
        ts_sim = current_time:dt:current_time + Nsim*dt;
        ts = [ts, ts_sim(2:end)];
        lengths = [lengths, length(ts)];
        current_time = current_time + Nsim*dt;
        x_traj = [x_traj, X_sim(:, 2:end)];
        u_traj = [u_traj, u_sim(:, 2:end)];
        P_traj = cat(3, P_traj, Pmat_sol(:, :, 2:end));
        ee_traj = [ee_traj, EE_ref_sol(:, 2:end)];
        P_EEPos = [P_EEPos, result.P_EEPos(:, 2:Nsim+1)];
        P_EEVel = [P_EEVel, result.P_EEVel(:, 2:Nsim+1)];

        % animate_trajectory(result);
        if abs(ee_traj(:, end) - [target_pos; 0; 0]) < [target_radius; target_radius; target_vel_accuracy; target_vel_accuracy]
            break
        end

        X_init = x_traj(:,end);
        P_init = P_0;
        u_init = u_traj(:,end);
        hot_start = true;
        max_iter = 1000;
        hs_X = result.X;
        hs_u = result.e_ff';
        hs_M = result.M;
        hs_T = result.time(end);
        dX_init = result.functions.f_forwardMusculoskeletalDynamics(X_init, u_init, 0);
    end
    if i == max_rollouts
        disp("Error: Max rollouts reached");
    end

    result.e_ff = u_traj';
    result.X = x_traj;
    result.q = x_traj(1:2, :)';
    result.qdot = x_traj(3:4, :)';
    result.Pmat = P_traj;
    result.time = ts;
    result.auxdata = result.auxdata;
    result.functions = result.functions;
    result.EEPos = ee_traj(1:2, :)';
    result.EEVel = ee_traj(3:4, :)';
    result.P_EEPos = P_EEPos;
    result.P_EEVel = P_EEVel;

end


% %% plot control inputs
% titles = {'m1 - Brachialis','m2 - Lateral triceps','m3 - anterior deltoid','m4 - posterior deltoid','m5 - biceps short','m6 - triceps long'};
% ymax = max(max(u_traj));
% figure;
% for i = 1:6
%     subplot(3,2,i)
%     hold on;
%     stairs(ts, u_traj(i,:), 'r', 'LineWidth', 2);
%     title(titles(i))
%     ylim([-0.1 0.1]);
%     xlim([0 ts(end)]);
%     xlabel('Time (s)');
%     ylabel('Activation');
%     % legend('Planned Trajectory (feedforward)', 'Actual Trajectory (feedback)', 'Actual $\pm$ std dev','Interpreter','latex');
% end

