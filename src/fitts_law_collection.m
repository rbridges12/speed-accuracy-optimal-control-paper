%% script to collect data for Fitts' Law

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
target_vel_accuracy = 0.2; % 95% confidence interval for final velocity radius
k_u = 1; % control effort weight
k_t = 10; % duration weight
% filename = "fitts_law_trials.mat";
filename = "fitts_law_trials_dfw.mat";

% cell array to store each data point
% each row is: [q_init (1-2), target_pos (3-4), target_radius (5), distance (6), time (7), failures (8), max_vel (9), max_vel_time (10)]
try
    load(filename);
catch
    trials = {};
end
% initial_pos = [0.0; 0.3];
X_init = [0.4061; 2.1532; 0; 0]; % short
% X_init = [0.3; 1.0; 0; 0]; % long
target_pos = [-0.1; .45]; % short
% target_pos = [-0.21; .50]; % long
% target_radius = 0.03; % 95% confidence interval for final position radius
target_radii = [0.015, 0.02, 0.03, 0.05, 0.06, 0.07, 0.08, 0.1];
q_inits = [[0.4061; 2.1532], [0.35; 2], [0.3; 1.5], [0.3; 2.5]];
colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'w'];
% times = [];w
% distances = [];
% radii = [];
% legend_labels = {};
num_trials = 5;
%%

figure
tiledlayout(2,1);
nexttile;
title('Normalized Velocity of End Effector')
xlabel('Normalized Time');
ylabel('Normalized Velocity');
hold on; grid on;
for i = 1:length(target_radii)
    for k = 1:length(q_inits)
        target_radius = target_radii(i);
        X_init = [q_inits(:,k); 0; 0];
        if trial_already_run(trials, X_init, target_pos, target_radius)
            continue
        end
        color = colors(i);
        max_vels = [];
        max_vel_times = [];
        trial_times = [];
        failures = 0;
        for j = 1:num_trials
            msg = sprintf("Running trial: i=%d, k=%d, j=%d", i, k, j);
            disp(msg);
            try
                result = nonlinear_mpc(N, Tsim, motor_noise_stddev, target_radius, target_vel_accuracy, k_u, k_t, X_init, target_pos);
                EE_vel = result.EEVel;
                norm_vel = vecnorm(EE_vel,2,2);
                [max_vel, max_vel_i] = max(norm_vel);
                max_vel_times = [max_vel_times, result.time(max_vel_i) / max(result.time)];
                max_vels = [max_vels, max_vel];
                trial_times = [trial_times, result.time(end)];

                normalized_vel = norm_vel./max(norm_vel);
                normalized_time = result.time./max(result.time);
                plot(normalized_time, normalized_vel, color, 'LineWidth', 2)
            catch E
                disp(E)
                failures = failures + 1;
            end
        end
        % legend_labels = [legend_labels, sprintf("radius: %f", target_radius)];
        movement_distance = norm(target_pos - EndEffectorPos(X_init(1:2), result.auxdata));
        avg_max_vel = mean(max_vels);
        avg_max_vel_time = mean(max_vel_times);
        trial_data = [X_init(1:2)', target_pos', target_radius, movement_distance, mean(trial_times), failures, avg_max_vel, avg_max_vel_time];
        trials = [trials; trial_data];
        save(filename, "trials");
        % distances = [distances, movement_distance];
        % radii = [radii, target_radius];
        % times = [times, mean(trial_times)];
    end
end
hold off

%% generate Fitts' Law plot

%  extract data from cell array
load(filename);
radii = [];
distances = [];
times = [];
for i = 1:length(trials)
    trial = trials{i};
    failures = trial(8);
    if failures >= num_trials - 1
        continue
    end
    radii = [radii, trial(5)];
    distances = [distances, trial(6)];
    times = [times, trial(7)];
end

% difficulty index = log base 2 of distance divided by target width
difficulties = log2(distances ./ (2 * radii));

% fit data to Fitts' Law: reaching time = a + b * difficulty index
A = [ones(length(difficulties), 1), difficulties'];
b = times';
% b = ones(8,1);
x = A\b;
a = x(1);
b = x(2);
fit_times = a + b * difficulties;

% nexttile;
figure;
title("Fitts' Law Fit, a = " + a + ", b = " + b);
xlabel('Difficulty Index');
ylabel('Movement Duration (s)');
hold on; grid on;
plot(difficulties, times, 'bo');
plot(difficulties, fit_times, 'r', 'LineWidth', 2);
legend('Data', 'Fitts'' Law Fit', 'location', 'northwest');

function result = trial_already_run(trials, X_init, target_pos, target_radius)
    for i = 1:length(trials)
        trial = trials{i};
        if isequal(trial(1:4), [X_init', target_pos', target_radius])
            result = true;
            return
        end
    end
    result = false;
end