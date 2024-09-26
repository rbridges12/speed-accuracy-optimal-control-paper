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
X_init = [0.2; 2; 0; 0];
% X_init = [ik_opt([0;0]); 0; 0];
P_0 = diag([1e-4; 1e-4; 1e-7; 1e-7]);
target_pos = [0.1; .6];
% target_pos = [0;.3];
target_radius = 0.05; % 95% confidence interval for final position radius
target_vel_accuracy = 0.2; % 95% confidence interval for final velocity radius
k_u = 1; % control effort weight
k_t = 1; % duration weight
color = 'b';

%% velocity profiling
max_vels = [];
max_vel_times = [];
times = [];
figure
hold on
xlabel('Normalized Time')
ylabel('Normalized Velocity')
for i = 1:5
    %result = nonlinear_mpc(N, Tsim, motor_noise_stddev, 0.1, target_vel_accuracy, 1, 1, X_init, target_pos);
    result_fast = nonlinear_mpc(N, Tsim, motor_noise_stddev, 0.035, target_vel_accuracy, 1, 100, [3.1416; 3.1416; 0; 0], [0.1; 0.35]);
    EE_vel = result.EEVel;
    norm_vel = vecnorm(EE_vel,2,2);
    [max_vel, max_vel_i] = max(norm_vel);
    max_vel_times = [max_vel_times, result.time(max_vel_i) / max(result.time)];
    max_vels = [max_vels, max_vel];
    times = [times, result.time(end)];

    normalized_vel = norm_vel./max(norm_vel);
    normalized_time = result.time./max(result.time);

    plot(normalized_time, normalized_vel, color, 'LineWidth', 2)
    % plot(normalized_time, norm_vel, 'LineWidth', 2)
    movement_distance = norm(target_pos - EndEffectorPos(X_init(1:2), result.auxdata));
    movement_time = mean(times);
    title('Normalized Velocity of End Effector')
    % t = sprintf("ku: %f, kt: %f, radius: %f, distance: %f, time: %f", k_u, k_t, target_radius, movement_distance, movement_time);
    % title(t)
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
N = 40; % number of discretized nodes
motor_noise_stddev = 0.036; % motor noise standard deviation
X_init = [0.2; 2; 0; 0];
% X_init = [ik_opt([0;0]); 0; 0];
P_init = diag([1e-4; 1e-4; 1e-7; 1e-7]);
target_pos = [0.1; .6];
% target_pos = [0;.3];
target_radius = 0.05; % 95% confidence interval for final position radius
target_vel_accuracy = 0.2; % 95% confidence interval for final velocity radius
k_u = 1; % control effort weight
k_t = 1; % duration weight

% ee_init_x = [0,0]
% ee_init_y = [0,.3]-[.1,.6]
% k_u = .05, k_t = 1
% r = .3, .4
%slow {[0.2 2 0.1000 0.6000 0.0500 ku=1 kt=1]}
%fast [     0.3500 2 0.1000 0.5000 0.1000 0.1699 0.2521 0 1.0458 0.6750 ku=1 kt=20]}

result = optimization_6muscles(N, motor_noise_stddev, target_radius, target_vel_accuracy, k_u, k_t, X_init, 0, 0, P_init, target_pos, false, [], [], [], 0, 3000);

%%
    result_fast = optimization_6muscles(N, motor_noise_stddev, 0.1, target_vel_accuracy, 1, 20, [0.35; 2; 0; 0], 0, 0, P_init, [0.1; 0.5], false, [], [], [], 0, 3000);
%%
    EE_vel = result.EEVel;
    norm_vel = vecnorm(EE_vel,2,2);
    [max_vel, max_vel_i] = max(norm_vel);
    max_vel_time = result.time(max_vel_i) / max(result.time);

    normalized_vel = norm_vel./max(norm_vel);
    normalized_time = result.time./max(result.time);

    EE_vel = result_fast.EEVel;
    norm_vel = vecnorm(EE_vel,2,2);
    [fast_max_vel, fast_max_vel_i] = max(norm_vel);
    max_vel_time_fast = result_fast.time(fast_max_vel_i) / max(result_fast.time);

    fast_normalized_vel = norm_vel./max(norm_vel);
    fast_normalized_time = result_fast.time./max(result_fast.time);

    figure
    blue = [0 0.4470 0.7410];
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'defaultLegendInterpreter','latex');
    alpha = 0.7;
    Darkgrey = [.25 .25 .25];
    DarkBlue = [0 .2 .4];
    DarkBlue = [DarkBlue, alpha];
    Azure = [53, 112, 188]/255;
    Azure = [Azure, alpha];
    VermillionRed = [156,31,46]/255;
    VermillionRed = [VermillionRed, alpha];
    DupontGray = [144,131,118]/255;
    grid on; hold on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14, "XMinorTick", "on", "YMinorTick", "on");
    plot(normalized_time, normalized_vel, 'Color', blue, 'LineWidth', 2)
    plot([max_vel_time max_vel_time], [0 1], 'Color', blue, 'LineStyle', '--', 'LineWidth', 1)
    plot(fast_normalized_time, fast_normalized_vel, 'Color', VermillionRed, 'LineWidth', 2)
    plot([max_vel_time_fast max_vel_time_fast], [0 1], 'Color', VermillionRed, 'LineStyle', '--', 'LineWidth', 1)
    
    % movement_distance = norm(target_pos - EndEffectorPos(X_init(1:2), result.auxdata));
    % movement_time = mean(times);
    % title('Normalized Velocity of End Effector')
    xlabel('Normalized Time', 'Interpreter', 'latex');
    ylabel('Normalized Velocity', 'Interpreter', 'latex');
    legend("slow", "", "fast", "", 'location', 'northwest');
    max_vel_time
    max_vel
    max_vel_time_fast
    fast_max_vel
    result.time(end)
    result_fast.time(end)
    slow_dist = norm(target_pos - EndEffectorPos(X_init(1:2), result.auxdata))
    fast_dist = norm([0.1; 0.5] - EndEffectorPos([0.35; 2], result.auxdata))
    % title(t)

%%
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