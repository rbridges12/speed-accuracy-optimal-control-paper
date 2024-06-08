auxdata = initializeModelParameters();
auxdata.N = N; % number of discretized nodes
nStates = 4; auxdata.nStates = nStates; % [q1, q2, qdot1, qdot2]
nControls = auxdata.nMuscles;

x = zeros(nStates,1);
P = diag([1e-2; 1e-2; 1e-3; 1e-3]);
kappa = 2;
% f = @(x) [x(1) * cos(x(2)); x(1) * sin(x(2))];
dt = 0.3;
u = ones(nControls,1)*0.5;
f = @(x) x + forwardMusculoskeletalDynamics_motorNoise(x, u, 0, zeros(auxdata.nMotorNoises,1), auxdata) * dt;


% sigma points around the reference point
L = sqrt(nStates + kappa) * chol(P, 'lower');
Y = x(:, ones(1, numel(x)));
X = [x, Y + L, Y - L];
w = zeros(2 * nStates + 1, 1);
w(1) = kappa / (nStates + kappa);
w(2:end) = 1 / (2*(nStates + kappa));

% % compute sample mean and covariance
mean = 0;
Y = [];
for i = 1:2*nStates+1
    Y(:,i) = f(X(:,i));
    mean = mean + w(i) * Y(:,i);
end
Cov = (Y - mean) * diag(w) * (Y - mean)';
Cov_xy = (X - x) * diag(w) * (Y - mean)';

%% arm visualization
figure
subplot(1,2,1)
hold on;
for i = 1:size(X, 2)
    draw_arm(X(1,i), X(2,i), auxdata.l1, auxdata.l2, 2, 2, 'b-', 'co', 'ro', false)
end
axis equal 
axis([-0.5 0.5 -0.1 0.8])
xlabel("X Position (m)")
ylabel("Y Position (m)")

subplot(1,2,2)
hold on;
for i = 1:size(Y, 2)
    draw_arm(Y(1,i), Y(2,i), auxdata.l1, auxdata.l2, 2, 2, 'b-', 'co', 'ro', false)
end
axis equal 
axis([-0.5 0.5 -0.1 0.8])
xlabel("X Position (m)")
ylabel("Y Position (m)")

%% 2D visualization
fsize = 22; % font size
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
green = [0.2980 .6 0];
% crimson = [220,20,60]/255; 
% darkblue = [0 .2 .4];
Darkgrey = [.25 .25 .25];
% darkgrey = [.35 .35 .35];
% lightgrey = [.7 .7 .7];
% Lightgrey = [.9 .9 .9];
VermillionRed = [156,31,46]/255;
DupontGray = [144,131,118]/255;
% Azure = [53, 112, 188]/255;
% purple = [178, 102, 255]/255;
% orange = [255,110,0]/255;

% create confidence ellipse
% first create points from a unit circle
phi = (-pi:.01:pi)';
circle = [cos(phi), sin(phi)];
% Chi-squared 2-DOF 95% percent confidence (0.05): 5.991
scale = sqrt(5.991);
% apply the transformation and scale of the covariance
ellipse_polar = (scale * chol(P,'lower') * circle' + x)';
ellipse_cartesian = (scale * chol(Cov,'lower') * circle' + mean)';

% generate samples for both polar and cartesian coordinates
s_polar = (chol(P,'lower') * randn(2,1000) + x)';
s_cartesian = zeros(size(s_polar));
for i = 1:size(s_polar,1)
    s_cartesian(i,:) = f(s_polar(i,:));
end

% plot in polar coordinates
figure; hold on; grid on
h = []; % plot handle
h{1} = plot(s_polar(:,1), s_polar(:,2), '.', 'color', DupontGray);
h{2} = plot(x(1), x(2), 'o', 'color', VermillionRed, 'markersize', 18);
h{3} = plot(ellipse_polar(:,1), ellipse_polar(:,2), 'color', VermillionRed, 'linewidth', 3);
h{4} = plot(X(1,:), X(2,:), '.', 'color', Darkgrey, 'markersize', 32);
xlabel('$r$', 'Interpreter','latex'); 
ylabel('$\theta$', 'Interpreter','latex');
legend([h{1}, h{2}, h{3}, h{4}], 'Samples', 'Mean', '$95\%$ Confidence Ellipse', 'Sigma Points', 'location', 'north outside')
text(1.75, 1.6, '$\kappa = 2$', 'fontsize',fsize, 'Interpreter','latex')
axis equal auto
set(gca,'fontsize',fsize)
set(gca,'TickLabelInterpreter','latex')
% figuresize(21,21,'cm')
% print -opengl -dpng -r600 ut_example_polar.png

% plot in Cartesian coordinates
figure; hold on; grid on
h = []; % plot handle
h{1} = plot(s_cartesian(:,1), s_cartesian(:,2), '.', 'color', DupontGray);
h{2} = plot(mean(1), mean(2), 'o', 'color', VermillionRed, 'markersize', 18);
h{3} = plot(ellipse_cartesian(:,1), ellipse_cartesian(:,2), 'color', VermillionRed, 'linewidth', 3);
h{4} = plot(Y(1,:), Y(2,:), '.', 'color', Darkgrey, 'markersize', 32);
xlabel('$x=r\cos(\theta)$', 'Interpreter','latex'); 
ylabel('$y=r\sin(\theta)$', 'Interpreter','latex');
legend([h{1}, h{2}, h{3}, h{4}], 'Samples', 'Mean', '$95\%$ Confidence Ellipse', 'Sigma Points', 'location', 'north outside')
text(1.6, 1.8, '$\kappa = 2$', 'fontsize',fsize, 'Interpreter','latex')
axis equal auto
set(gca,'fontsize',fsize)
set(gca,'TickLabelInterpreter','latex')
% figuresize(21,21,'cm')