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
N = 80; % number of discretized nodes
motor_noise_stddev = 0.03; % motor noise standard deviation
initial_pos = [0; 0.3];
target_pos = [-0.1; .45];
target_pos_accuracy = 0.03; % 95% confidence interval for final position radius
target_vel_accuracy = 0.1; % 95% confidence interval for final velocity radius
k_u = 0.1; % control effort weight
k_t = 1; % duration weight

result = optimization_6muscles(N, motor_noise_stddev, target_pos_accuracy, target_vel_accuracy, k_u, k_t, initial_pos, target_pos);

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
sigma_test(result);

