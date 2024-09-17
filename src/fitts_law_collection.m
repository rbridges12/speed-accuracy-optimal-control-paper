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

filename = "fitts_law_ethan_real_single_weight1.mat";
type = "single";
% filename = "fitts_law_ethans_real.mat";
% type = "mpc";

N = 40; % number of discretized nodes
Tsim = 0.1;
motor_noise_stddev = 0.036; % motor noise standard deviation
target_vel_accuracy = 0.2; % 95% confidence interval for final velocity radius
k_u = 1; % control effort weight
k_t = 5; % duration weight
color = "r";
num_trials = 1;
if type == "mpc"
    color = "b";
    num_trials = 3;
end

% cell array to store each data point
% each row is: [q_init (1-2), target_pos (3-4), target_radius (5), distance (6), time (7), failures (8), max_vel (9), max_vel_time (10)]
try
    load(filename);
catch
    trials = {};
end
% X_init = [0.4061; 2.1532; 0; 0]; % short
% X_init = [ik_opt([0; 0]); 0; 0];
% X_init = [0.3; 1.0; 0; 0]; % long
% target_pos = [-0.1; .45]; % short
% target_pos = [0; 0.3]; % short
% target_pos = [0.1; 0.3]; % short
% target_pos = [-0.21; .50]; % long
% target_radius = 0.03; % 95% confidence interval for final position radius
P_init = diag([1e-4; 1e-4; 1e-7; 1e-7]);
target_radii = [0.03, 0.04];
q_inits = [ik_opt([0; .3])];
target_ps = [[0; 0.35], [0; 0.36], [0; 0.38], [0; 0.4], [0; 0.45], [0; 0.5], [0; 0.55], [0; 0.6], [0.1; 0.35], [0.1; 0.36], [0.1; 0.38], [0.1; 0.4], [0.1; 0.45], [0.1; 0.5], [0.1; 0.55], [0.1; 0.6]];

%%
figure
title('Normalized Velocity of End Effector')
xlabel('Normalized Time');
ylabel('Normalized Velocity');
hold on; grid on;
for i = 1:length(target_radii)
    for k = 1:size(q_inits, 2)
        for l = 1:size(target_ps, 2)
            target_radius = target_radii(i);
            X_init = [q_inits(:,k); 0; 0];
            target_pos = target_ps(:,l);
            if trial_already_run(trials, X_init, target_pos, target_radius)
                continue
            end
            max_vels = [];
            max_vel_times = [];
            trial_times = [];
            failures = 0;
            for j = 1:num_trials
                msg = sprintf("Running trial: i=%d, k=%d, l=%d, j=%d", i, k, l, j);
                disp(msg);
                try
                    if type == "single"
                        result = optimization_6muscles(N, motor_noise_stddev, target_radius, target_vel_accuracy, k_u, k_t, X_init, zeros(6,1), zeros(4,1), P_init, target_pos, false, [], [], [], 0, 3000);
                    elseif type == "mpc"
                        result = nonlinear_mpc(N, Tsim, motor_noise_stddev, target_radius, target_vel_accuracy, k_u, k_t, X_init, target_pos);
                    end
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
            movement_distance = norm(target_pos - EndEffectorPos(X_init(1:2), result.auxdata));
            avg_max_vel = mean(max_vels);
            avg_max_vel_time = mean(max_vel_times);
            trial_data = [X_init(1:2)', target_pos', target_radius, movement_distance, mean(trial_times), failures, avg_max_vel, avg_max_vel_time];
            trials = [trials; trial_data];
            save(filename, "trials");
        end
    end
end
hold off

%% generate Fitts' Law plot
%  extract data from cell array
load(filename);
%load("fitts_law_trials_dfw.mat")
radii = [];
distances = [];
times = [];
for i = 1:length(trials)
    trial = trials{i};
    failures = trial(8);
    difficulty = log2(trial(6) ./ (2 * trial(5)));
    time = trial(7);
    if failures >= num_trials || difficulty < 0 || time > 1
        out = sprintf("skipped trial: failures=%d, difficulty=%f, time=%f, radius=%f, q_init=[%f, %f], target_pos=[%f, %f]", failures, difficulty, time, trial(5), trial(1), trial(2), trial(3), trial(4));
        disp(out);
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
x = A\b;
a = x(1);
b = x(2);
fit_times = a + b * difficulties;
mdl = fitlm(difficulties, times);
rsq = mdl.Rsquared.Ordinary;

figure;
title("Fitts' Law Fit, a = " + a + ", b = " + b + ", R^2 = " + rsq);
xlabel('Difficulty Index');
ylabel('Movement Duration (s)');
hold on; grid on;
plot(difficulties, times, 'bo');
plot(difficulties, fit_times, 'r', 'LineWidth', 2);
legend('Data', 'Fitts'' Law Fit', 'location', 'northwest');

function result = trial_already_run(trials, X_init, target_pos, target_radius)
    for i = 1:length(trials)
        trial = trials{i};
        if isequal(trial(1:5), [X_init(1:2)', target_pos', target_radius])
            result = true;
            return
        end
    end
    result = false;
end