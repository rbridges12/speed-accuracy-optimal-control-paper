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

N = 40; % number of discretized nodes
Tsim = 0.1;
motor_noise_stddev = 0.036; % motor noise standard deviation
% initial_pos = [0.0; 0.3];
% X_init = [0.4061; 2.1532; 0; 0]; % short
X_init = [0.3; 1.0; 0; 0]; % long
% target_pos = [-0.1; .45]; % short
target_pos = [-0.21; .50]; % long
target_radius = 0.03; % 95% confidence interval for final position radius
target_vel_accuracy = 0.2; % 95% confidence interval for final velocity radius
k_u = 1; % control effort weight
k_t = 5; % duration weight

%% velocity profiling
max_vels = [];
max_vel_times = [];
times = [];
figure
hold on
xlabel('Normalized Time')
ylabel('Normalized Velocity')
for i = 1:5
    result = nonlinear_mpc(N, Tsim, motor_noise_stddev, target_radius, target_vel_accuracy, k_u, k_t, X_init, target_pos);
    EE_vel = result.EEVel;
    norm_vel = vecnorm(EE_vel,2,2);
    [max_vel, max_vel_i] = max(norm_vel);
    max_vel_times = [max_vel_times, result.time(max_vel_i) / max(result.time)];
    max_vels = [max_vels, max_vel];
    times = [times, result.time(end)];

    normalized_vel = norm_vel./max(norm_vel);
    normalized_time = result.time./max(result.time);

    plot(normalized_time, normalized_vel, 'LineWidth', 2)
    % plot(normalized_time, norm_vel, 'LineWidth', 2)
    movement_distance = norm(target_pos - EndEffectorPos(X_init(1:2), result.auxdata));
    movement_time = mean(times);
    % title('Normalized Velocity of End Effector')
    t = sprintf("ku: %f, kt: %f, radius: %f, distance: %f, time: %f", k_u, k_t, target_radius, movement_distance, movement_time);
    title(t)
end
hold off
avg_max_vel = mean(max_vels)
avg_max_vel_time = mean(max_vel_times)

%% single MPC run
result = nonlinear_mpc(N, Tsim, motor_noise_stddev, target_radius, target_vel_accuracy, k_u, k_t, X_init, target_pos);

%% plot MPC trajectory
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

%% feedforward trajectory optimization
% set model and optimization parameters
N = 80; % number of discretized nodes
motor_noise_stddev = 0.036; % motor noise standard deviation
% initial_pos = [-0.2; 0.5];
X_init = [0.4061; 2.1532; 0; 0];
P_init = diag([1e-4; 1e-4; 1e-7; 1e-7]);
% initial_pos = q_init;
target_pos = [-0.1; .45];
target_radius = 0.05; % 95% confidence interval for final position radius
target_vel_accuracy = 0.2; % 95% confidence interval for final velocity radius
k_u = 0.1; % control effort weight
k_t = 0.1; % duration weight

result = optimization_6muscles(N, motor_noise_stddev, target_radius, target_vel_accuracy, k_u, k_t, X_init, P_init, target_pos);
animate_trajectory(result);

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

%% save ff data
fname = sprintf("ff_%f_%f_%f_%f_%s.mat", target_radius, target_vel_accuracy, k_u, k_t, datestr(now, 'dd-HH-MM'));
save(fname, 'result');

%% save MPC data
fname = sprintf("mpc_%f_%f_%f_%f_%s.mat", target_radius, target_vel_accuracy, k_u, k_t, datestr(now, 'dd-HH-MM'));
save(fname, 'result');