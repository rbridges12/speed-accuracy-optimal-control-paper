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

%% discretization test
% u = @(t) ones(6, 1) * 0.1;
u = @(t) interp1(result.time, result.e_ff, t)';
u_gt = u;

% simulate the continuous time dynamics
continuous_ode = @(t, x) forwardMusculoskeletalDynamics_motorNoise(x, u(t), 0, zeros(6, 1), result.auxdata);
[ts_c, xs_c] = ode45(continuous_ode, [0, result.time(end)], result.X(:, 1));

% simulate the discrete time dynamics
dt = result.time(end) / size(result.X, 2);
discrete_ode = @(x, u) discrete_dynamics(x, u, zeros(6, 1), result.auxdata, dt);
ts_d = 0:dt:result.time(end);
xs_d = zeros(4, length(ts_d));
xs_d(:, 1) = result.X(:, 1);
for i = 2:length(ts_d)
    xs_d(:, i) = discrete_ode(xs_d(:, i - 1), u(ts_d(i)));
end

% plot the results
figure;
subplot(2, 2, 1);
plot(ts_c, xs_c(:, 1), 'r', ts_d, xs_d(1, :), 'b--', 'LineWidth', 2);
title('q1');
legend('continuous', 'discrete');
grid on;
subplot(2, 2, 2);
plot(ts_c, xs_c(:, 2), 'r', ts_d, xs_d(2, :), 'b--', 'LineWidth', 2);
title('q2');
legend('continuous', 'discrete');
grid on;
subplot(2, 2, 3);
plot(ts_c, xs_c(:, 3), 'r', ts_d, xs_d(3, :), 'b--', 'LineWidth', 2);
title('qdot1');
legend('continuous', 'discrete');
grid on;
subplot(2, 2, 4);
plot(ts_c, xs_c(:, 4), 'r', ts_d, xs_d(4, :), 'b--', 'LineWidth', 2);
title('qdot2');
legend('continuous', 'discrete');
grid on;

%% linearization test
dt = 0.005; % discretization time step
lin_interval = 1; % time between updates of the linearization point

% simulate dynamics
f = @(x, u) discrete_dynamics(x, u, zeros(6, 1), result.auxdata, dt);
ts = 0:dt:result.time(end);
xs_d = zeros(4, length(ts));
xs_d(:, 1) = result.X(:, 1);
xs_dl = zeros(4, length(ts));
xs_dl(:, 1) = result.X(:, 1);
xstar = xs_dl(:, 1);
ustar = u(0);
for i = 2:length(ts)
    % nonlinear system
    xs_d(:, i) = f(xs_d(:, i - 1), u(ts(i)));
    
    % linearized system
    if mod(i, lin_interval / dt) == 0
        xstar = xs_dl(:, i - 1);
        ustar = u(ts(i));
    end
    [A, B] = finite_diff_jacobians(f, xstar, ustar);
    xs_dl(:, i) = f(xstar, ustar) + A * (xs_dl(:, i - 1) - xstar) + B * (u(ts(i)) - ustar);
end

% plot the results
figure;
subplot(2, 2, 1);
plot(ts, xs_d(1, :), 'b', ts, xs_dl(1, :), 'r--', 'LineWidth', 2);
title('q1');
legend('nonlinear', 'linearized');
grid on;
subplot(2, 2, 2);
plot(ts, xs_d(2, :), 'b', ts, xs_dl(2, :), 'r--', 'LineWidth', 2);
title('q2');
grid on;
subplot(2, 2, 3);
plot(ts, xs_d(3, :), 'b', ts, xs_dl(3, :), 'r--', 'LineWidth', 2);
title('qdot1');
grid on;
subplot(2, 2, 4);
plot(ts, xs_d(4, :), 'b', ts, xs_dl(4, :), 'r--', 'LineWidth', 2);
title('qdot2');
grid on;

%% QP formulation
dt = 0.005; % discretization time step
n = 4; l = 6; % state and control dimensions
N_h = 20; % prediction horizon
Q = eye(n); % state cost
R = 0.0001 * eye(l); % input cost
Qf = 100 * eye(n); % terminal state cost
x_init = result.X(:, 1); % initial state
x_target = result.X(:, end); % target state
u_min = 0.001; u_max = 1; % control bounds
q_min = 0; q_max = 180; % state bounds
  
% define block cost matrices
HQ = blkdiag(kron(eye(N_h - 1), Q), Qf);
HR = kron(eye(N_h - 1), R);

% linearized dynamics
f = @(x, u) discrete_dynamics(x, u, zeros(6, 1), result.auxdata, dt);
xstar = zeros(n, 1);
ustar = zeros(l, 1);
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
        x(:, 1) == x_init;

cvx_end

% plot the results
figure;
subplot(2, 2, 1);
plot(1:N_h, x(1, :), 'b', 'LineWidth', 2);
title('q1');
grid on;
subplot(2, 2, 2);
plot(1:N_h, x(2, :), 'b', 'LineWidth', 2);
title('q2');
grid on;
subplot(2, 2, 3);
plot(1:N_h, x(3, :), 'b', 'LineWidth', 2);
title('qdot1');
grid on;
subplot(2, 2, 4);
plot(1:N_h, x(4, :), 'b', 'LineWidth', 2);
title('qdot2');
grid on;

figure;
hold on;
for i = 1:l
    plot(1:N_h - 1, u(i, :), 'LineWidth', 2);
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
