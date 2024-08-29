%% Main script for the optimization of the arm model

% include functions in subdirectories
addpath("~/casadi-3.6.5")
addpath("./forwardSim")
addpath("./Muscle_LMT_dM")
addpath("./MuscleModel")
addpath("./ArmModel")
addpath("./MusculoskeletalDynamics")
addpath("./Integrator")
addpath("./plotFunctions")

% set model and optimization parameters
N = 40; % number of discretized nodes
motor_noise_stddev = 0.036; % motor noise standard deviation
% initial_pos = [0.0; 0.3];
X_init = [0.4061; 2.1532; 0; 0];
P_init = diag([1e-4; 1e-4; 1e-7; 1e-7]);
target_pos = [-0.1; .45];
target_radius = 0.1; % 95% confidence interval for final position radius
target_vel_accuracy = 0.2; % 95% confidence interval for final velocity radius
k_u = 0.1; % control effort weight
k_t = 0.1; % duration weight
Tsim = 0.05;

current_time = 0;
dX_init = zeros(4, 1);
u_init = zeros(6, 1);
hot_start = false;
ts = [];
ts_plans = [];
X_plans = [];
lengths = [];
x_traj = [];
u_traj = [];
P_traj = [];
ee_traj = [];

for i = 1:50
    result = optimization_6muscles(N, motor_noise_stddev, target_radius, target_vel_accuracy, k_u, k_t, X_init, u_init, hot_start, dX_init, P_init, target_pos);

    dt = result.time(2) - result.time(1);
    Nsim = min(ceil(Tsim / dt), N);
    u_sim = result.e_ff';
    u_sim = u_sim(:, 1:Nsim+1);
    [X_sim, ~, EE_ref_sol, Pmat_sol] = forwardSim_ode(result.X(:,1), result.Pmat(:,:,1) ,result.e_ff', Nsim*dt, Nsim, result.auxdata, result.functions);

    ts_planned = result.time + current_time;
    ts_plans = [ts_plans; ts_planned];
    X_plans = cat(3, X_plans, result.X);
    ts_sim = current_time:dt:current_time + Nsim*dt;
    ts = [ts, ts_sim];
    lengths = [lengths, length(ts)];
    current_time = current_time + Nsim*dt;
    x_traj = [x_traj, X_sim];
    u_traj = [u_traj, u_sim];
    P_traj = cat(3, P_traj, Pmat_sol);
    ee_traj = [ee_traj, EE_ref_sol];

    % animate_trajectory(result);

    if ee_traj(1:2, end) - target_pos < 0.0001
        break
    end
    X_init = x_traj(:,end);
    P_init = P_traj(:,:,end);
    u_init = u_traj(:,end);
    hot_start = true;
    dX_init = result.functions.f_forwardMusculoskeletalDynamics(X_init, u_init, 0);
end

%% plot trajectory
for i = 1:length(lengths)
    titles = {'q1','q2','qdot1','qdot2'};
    ylabels = {"Shoulder Angle (rad)", "Elbow Angle (rad)", "Shoulder Velocity (rad/s)", "Elbow Velocity (rad/s)"};
    tiledlayout(2, 2);
    title("MPC Trajectory");
    for j = 1:4
        nexttile; cla;
        hold on; grid on;
        % plot(result.time, result.X(j, :), 'r', 'LineWidth', 2);
        % plot(ts, x_traj(j, :), 'b', 'LineWidth', 2);
        plot(ts(1:lengths(i)), x_traj(j, 1:lengths(i)), 'b', 'LineWidth', 2);
        plot(ts_plans(i, :), squeeze(X_plans(j, :, i)), 'g--', 'LineWidth', 2);
        title(titles(j));
        xlabel('Time (s)');
        ylabel(ylabels{j});
    end
    % legend("Reference", "Actual Trajectory", "Planned Trajectory", "location", "west outside");
    legend("Actual Trajectory", "Planned Trajectory", "location", "west outside");
    drawnow
    pause(1);
end

%% plot control inputs
titles = {'m1 - Brachialis','m2 - Lateral triceps','m3 - anterior deltoid','m4 - posterior deltoid','m5 - biceps short','m6 - triceps long'};
ymax = max(max(u_traj));
figure;
for i = 1:6
    subplot(3,2,i)
    hold on;
    stairs(ts, u_traj(i,:), 'r', 'LineWidth', 2);
    title(titles(i))
    ylim([-0.1 0.1]);
    xlim([0 ts(end)]);
    xlabel('Time (s)');
    ylabel('Activation');
    % legend('Planned Trajectory (feedforward)', 'Actual Trajectory (feedback)', 'Actual $\pm$ std dev','Interpreter','latex');
end


%% print covariances
result.final_cost;
2.*sqrt(result.final_pos_cov);
2.*sqrt(result.final_vel_cov);

covs = result.P_EEPos(:, end);
P_EE_final = [covs(1) covs(2); covs(2) covs(3)]
target_P = [(target_pos_accuracy/2)^2 0; 0 (target_pos_accuracy/2)^2]

%%
animate_trajectory(result);

%%
plotNormalizedVelocity(result);

%%
plotResults(result);

%% print trajectory duration
format long
duration = result.time(end)

%%
distributions_test(result);

%% save data
save('result.mat', 'result');