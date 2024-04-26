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
wM_std = 0.01; % motor noise standard deviation
EE_init = [0; 0.3];
EE_target = [0; .45];
final_pos_variance_95 = 0.03 ; % 95% confidence interval for final position radius
final_vel_variance_95 = 0.1; % 95% confidence interval for final velocity radius
k_u = .05; % control effort weight
k_t = 1; % duration weight

result = optimization_6muscles(N, wM_std, final_pos_variance_95, final_vel_variance_95, k_u, k_t, EE_init, EE_target);

%%
result.final_cost;
2.*sqrt(result.final_pos_cov);
2.*sqrt(result.final_vel_cov);

covs = result.P_EEPos(:, end);
P_EE_final = [covs(1) covs(2); covs(2) covs(3)]
target_P = [(final_pos_variance_95/2)^2 0; 0 (final_pos_variance_95/2)^2]

% figure
% hold on
% error_ellipse(P_EE_final, [0; 0], 0.95, "k:");
% error_ellipse(target_P, [0; 0], 0.95, "r:");
% legend("actual", "target")
% hold off

%%
animate_trajectory(result);

%%
plotNormalizedVelocity(result);

%%
plotResults(result);

%% 
format long
result.time(length(result.time))

