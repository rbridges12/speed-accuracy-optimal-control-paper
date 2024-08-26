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
motor_noise_stddev = 0.036; % motor noise standard deviation
% initial_pos = [0.0; 0.3];
initial_q = [0.4061; 2.1532];
target_pos = [-0.1; .45];
target_radius = 0.023; % 95% confidence interval for final position radius
target_vel_accuracy = 0.1; % 95% confidence interval for final velocity radius
k_u = 0.1; % control effort weight
k_t = 0; % duration weight
Tsim = 0.1;

x_traj = [];
u_traj = [];
P_traj = [];
ee_traj = [];

for i = 1:10
    result = optimization_6muscles(N, motor_noise_stddev, target_radius, target_vel_accuracy, k_u, k_t, initial_pos, target_pos);

    dt = result.time(2) - result.time(1);
    Nsim = ceil(Tsim / dt);
    u_sim = result.e_ff';
    u_sim = u_sim(1:Nsim+1);
    [X_sim, ~, EE_ref_sol, Pmat_sol] = forwardSim_ode(result.X(:,1), result.Pmat(:,:,1) ,result.e_ff', Nsim*dt, Nsim, result.auxdata, result.functions);

    x_traj = [x_traj, X_sim];
    u_traj = [u_traj, u_sim];
    P_traj = [P_traj, Pmat_sol];
    ee_traj = [ee_traj, EE_ref_sol];

    initial_pos = ee_traj(1:2, end);
    animate_trajectory(result);
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