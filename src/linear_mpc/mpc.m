% include functions in subdirectories
addpath("~/casadi-3.6.5")
addpath("./forwardSim")
addpath("./Muscle_LMT_dM")
addpath("./MuscleModel")
addpath("./ArmModel")
addpath("./MusculoskeletalDynamics")
addpath("./Integrator")
addpath("./plotFunctions")

load('result.mat');

dt = 0.005; % discretization time step
n = 4; l = 6; % state and control dimensions
t_h = 0.5; N_h = round(t_h / dt); % prediction horizon
u_steps = 1; % control steps per prediction horizon
iterations = 150;
Q = diag([10, 10, 1, 1]); % state cost
R = 100 * eye(l); % input cost
Qf = diag([1000, 1000, 100, 100]); % terminal state cost
x_init = result.X(:, 1); % initial state
x_target = result.X(:, end); % target state
u_min = 0.001; u_max = 1; % control bounds
q_min = 0; q_max = 180; % state bounds
  
f = @(x, u) discrete_dynamics(x, u, zeros(6, 1), result.auxdata, dt);

% define block cost matrices
HQ = blkdiag(kron(eye(N_h - 1), Q), Qf);
HR = kron(eye(N_h - 1), R);

% MPC loop
figure;
x_traj = [x_init];
x_lin_traj = [x_init];
u_traj = [zeros(l, 1)];
for i = 1:iterations
    % linearize dynamics around current state and control
    xstar = x_traj(:, end);
    ustar = u_traj(:, end);
    [A, B] = finite_diff_jacobians(f, xstar, ustar);

    % rollout horizon trajectory using convex optimization
    cvx_begin quiet
        variable x(n, N_h);
        variable u(l, N_h - 1);
        x_error = x - repmat(x_target, 1, N_h);
        minimize(quad_form(x_error(:), HQ) + quad_form(u(:), HR))

        subject to
            for i = 1:N_h - 1
                x(:, i + 1) == f(xstar, ustar) + A * (x(:, i) - xstar) + B * (u(:, i) - ustar);
                q_min <= x(1:2, i) <= q_max;
                u_min <= u(:, i) <= u_max;
            end
            q_min <= x(1:2, end) <= q_max;

            x(:, 1) == x_traj(:, end);
    cvx_end
    t_current = (size(x_traj, 2) - 1) * dt;
    h_ts = t_current:dt:t_current + (size(x, 2) - 1) * dt;

    % execute the first few steps of the trajectory on the "actual" system (nonlinear dynamics)
    x_lin_traj(:, end) = x_traj(:, end);
    for j = 1:u_steps
        x_next = f(x_traj(:, end), u(:, j));
        x_traj = [x_traj x_next];
        u_traj = [u_traj u(:, j)];

        x_lin_next = f(xstar, ustar) + A * (x_lin_traj(:, end) - xstar) + B * (u(:, j) - ustar);
        x_lin_traj = [x_lin_traj x_lin_next];
    end

    % live animation
    ts = 0:dt:(size(x_traj, 2) - 1) * dt;
    titles = {'q1','q2','qdot1','qdot2'};
    ylabels = {"Shoulder Angle (rad)", "Elbow Angle (rad)", "Shoulder Velocity (rad/s)", "Elbow Velocity (rad/s)"};
    tiledlayout(2, 2);
    title("MPC Trajectory");
    for j = 1:4
        nexttile; cla;
        hold on; grid on;
        plot(result.time, result.X(j, :), 'r', 'LineWidth', 2);
        plot(ts, x_traj(j, :), 'b', 'LineWidth', 2);
        plot(h_ts, x(j, :), 'g--', 'LineWidth', 2);
        title(titles(j));
        xlabel('Time (s)');
        ylabel(ylabels{j});
    end
    legend("Reference", "Actual Trajectory", "Planned Trajectory", "location", "west outside");
    drawnow
end

%% plot the results
ts = 0:dt:(size(x_traj, 2) - 1) * dt;
figure;
titles = {'q1','q2','qdot1','qdot2'};
for i = 1:4
    subplot(2, 2, i);
    hold on; grid on;
    plot(result.time, result.X(i, :), 'r', 'LineWidth', 2);
    plot(ts, x_traj(i, :), 'b', 'LineWidth', 2);
    plot(ts, x_lin_traj(i, :), 'g--', 'LineWidth', 2);
    title(titles(i));
    legend("Nonlinear TrajOpt", "MPC Trajectory", "MPC Horizon");
    xlabel('Time (s)');
end

figure;
titles = {'m1 - Brachialis','m2 - Lateral triceps','m3 - anterior deltoid','m4 - posterior deltoid','m5 - biceps short','m6 - triceps long'};
for i = 1:6
    subplot(3, 2, i);
    hold on;
    stairs(result.time, result.e_ff(:, i), 'r', 'LineWidth', 2);
    stairs(ts, u_traj(i, :), 'b', 'LineWidth', 2);
    title(titles(i));
    ylim([-0.1 0.1]);
    xlim([0 result.time(end)]);
    xlabel('Time (s)');
    ylabel('Activation');
    legend("Nonlinear TrajOpt", "Linearized MPC");
    grid on;
end

%%
function Xk1 = discrete_dynamics(Xk, uk, wMk, auxdata, dt)
    % use Runge-Kutta 4th order to discretize the continuous time dynamics
    f = @(x, u) forwardMusculoskeletalDynamics_motorNoise(x, u, 0, zeros(6, 1), auxdata);
    k1 = f(Xk, uk);
    k2 = f(Xk + dt / 2 * k1, uk);
    k3 = f(Xk + dt / 2 * k2, uk);
    k4 = f(Xk + dt * k3, uk);
    Xk1 = Xk + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
end

function [A, B] = finite_diff_jacobians(f, x, u)
    % compute the jacobians of the function f at x with respect to x and u
    n = length(x);
    l = length(u);
    A = zeros(n, n);
    B = zeros(n, l);
    eps = 1e-6;
    for i = 1:n
        x1 = x;
        x1(i) = x1(i) + eps;
        A(:, i) = (f(x1, u) - f(x, u)) / eps;
    end
    for i = 1:l
        u1 = u;
        u1(i) = u1(i) + eps;
        B(:, i) = (f(x, u1) - f(x, u)) / eps;
    end
end