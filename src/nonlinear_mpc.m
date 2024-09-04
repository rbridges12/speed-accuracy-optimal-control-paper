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
P_0 = diag([1e-4; 1e-4; 1e-7; 1e-7]);
target_pos = [-0.1; .45];
target_radius = 0.04; % 95% confidence interval for final position radius
target_vel_accuracy = 0.2; % 95% confidence interval for final velocity radius
k_u = 0.0; % control effort weight
k_t = 10; % duration weight
Tsim = 0.1;

P_init = P_0;
current_time = 0;
dX_init = zeros(4, 1);
u_init = zeros(6, 1);
hot_start = false;
hs_X = [];
hs_u = [];
hs_M = [];
hs_T = 0;
max_iter = 2000;
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

for i = 1:50
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

ts(end)

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
fname = sprintf("result_%f_%f_%f_%f_%s.mat", target_radius, target_vel_accuracy, k_u, k_t, datestr(now, 'dd-HH-MM'));
save(fname, 'result');