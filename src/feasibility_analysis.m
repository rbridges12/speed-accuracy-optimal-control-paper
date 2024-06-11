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
initial_pos = [0; 0.3];
target_pos = [-0.1; .45];
target_vel_accuracy = 0.1; % 95% confidence interval for final velocity radius
k_u = 0.1; % control effort weight
k_t = 1; % duration weight

motor_noise_stddev = linspace(0, 0.04, 10);
target_pos_accuracy = linspace(0.001, 0.1, 10);
feasibility_grid = zeros(length(motor_noise_stddev), length(target_pos_accuracy));
errors = [];

for i = 1:length(motor_noise_stddev)
    for j = 1:length(target_pos_accuracy)
        try
            tic;
            result = optimization_6muscles(N, motor_noise_stddev(i), target_pos_accuracy(j), target_vel_accuracy, k_u, k_t, initial_pos, target_pos);
            t_elapsed = toc;
            feasibility_grid(i, j) = t_elapsed;
        catch E
            disp(E)
            errors = [errors; E];
            feasibility_grid(i, j) = -1;
        end
    end
end
save("feasibility1.mat", "feasibility_grid", "motor_noise_stddev", "target_pos_accuracy", "errors")

%% Analysis
load("feasibility1.mat")
feasibility_grid(feasibility_grid == -1) = NaN;
xdata = compose("%.3f", target_pos_accuracy);
ydata = compose("%.3f", motor_noise_stddev);
heatmap(xdata, ydata, feasibility_grid, "MissingDataLabel", "Failed", "CellLabelFormat", "%.f", "GridVisible", "off")
% heatmap(target_pos_accuracy, motor_noise_stddev, feasibility_grid, "MissingDataLabel", "Optimization Failed", "MissingDataColor", [255, 180, 180]/255)
% colormap parula
xlabel("Target Size (m)")
ylabel("Motor Noise Standard Deviation")
title("Optimization Time (s) vs. Target Size and Motor Noise")