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
wM_std = 0.05; % motor noise standard deviation
EE_init = [0; 0.3];
EE_target = [0; 0.45];
final_pos_variance_95 = 0.1; % 95% confidence interval for final position
final_vel_variance_95 = 0.1; % 95% confidence interval for final velocity
k_u = 1; % control effort weight
k_t = 0; % duration weight

result = optimization_6muscles(N, wM_std, final_pos_variance_95, final_vel_variance_95, k_u, k_t, EE_init, EE_target);

result.final_cost;
2.*sqrt(result.final_pos_cov);
2.*sqrt(result.final_vel_cov);

%%
animate_trajectory(result);

%%
plotResults(result);
