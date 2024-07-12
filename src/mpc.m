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
N_h = 20; % prediction horizon
u_steps = 3; % control steps per prediction horizon
Q = eye(n); % state cost
R = 0.0001 * eye(l); % input cost
Qf = 100 * eye(n); % terminal state cost
x_init = result.X(:, 1); % initial state
x_target = result.X(:, end); % target state
u_min = 0.001; u_max = 1; % control bounds
q_min = 0; q_max = 180; % state bounds
x_traj = [x_init];
u_traj = [zeros(l, 1)];
  
f = @(x, u) discrete_dynamics(x, u, zeros(6, 1), result.auxdata, dt);

% define block cost matrices
HQ = blkdiag(kron(eye(N_h - 1), Q), Qf);
HR = kron(eye(N_h - 1), R);

for i = 1:100
    % linearized dynamics
    xstar = x_traj(:, end);
    ustar = u_traj(:, end);
    [A, B] = finite_diff_jacobians(f, xstar, ustar);

    cvx_begin
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
            x(:, 1) == x_traj(:, end);
    cvx_end

    for j = 1:u_steps
        x_next = f(x_traj(:, end), u(:, j));
        x_traj = [x_traj x_next];
        u_traj = [u_traj u(:, j)];
    end
end

%% plot the results
figure;
subplot(2, 2, 1);
plot(1:size(x_traj, 2), x_traj(1, :), 'b', 'LineWidth', 2);
title('q1');
grid on;
subplot(2, 2, 2);
plot(1:size(x_traj, 2), x_traj(2, :), 'b', 'LineWidth', 2);
title('q2');
grid on;
subplot(2, 2, 3);
plot(1:size(x_traj, 2), x_traj(3, :), 'b', 'LineWidth', 2);
title('qdot1');
grid on;
subplot(2, 2, 4);
plot(1:size(x_traj, 2), x_traj(4, :), 'b', 'LineWidth', 2);
title('qdot2');
grid on;

figure;
hold on;
for i = 1:l
    plot(1:size(u_traj, 2), u_traj(i, :), 'LineWidth', 2);
end
title('u');
grid on;


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