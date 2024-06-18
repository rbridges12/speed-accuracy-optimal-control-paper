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
u = @(t) ones(6, 1) * 0.1;

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


function Xk1 = discrete_dynamics(Xk, uk, wMk, auxdata, dt)
    % use Runge-Kutta 4th order to discretize the continuous time dynamics
    f = @(x, u) forwardMusculoskeletalDynamics_motorNoise(x, u, 0, zeros(6, 1), auxdata);
    k1 = f(Xk, uk);
    k2 = f(Xk + dt / 2 * k1, uk);
    k3 = f(Xk + dt / 2 * k2, uk);
    k4 = f(Xk + dt * k3, uk);
    Xk1 = Xk + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
end
