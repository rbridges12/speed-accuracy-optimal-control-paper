%  This file runs a comparison between the shooting and direct collocation
% implementation of the covariance matrix propagation, and this for a
% specific noise settings and target

% Note that the iterations for the shooting simulation run faster than for
% the DC simulation clearly for the BAR and CIRCLE case. However, for the
% OBSTACLE case the iteration times become more similar. This has to do
% with the increase in the number of path constraints in the OBSTACLE case
% where the sparsity is exploited in the DC implementation.

% The DC problems seem to take longer to converge (especially when the
% problems are initialized with a feasible and relatively good initial guess).
% The computational efficiency of the two implementations could probably be
% more equivalent if for the DC implementation some regularization is
% introduced. However, despite usually with regularization we can get equivalent or very similar results, 
% this would off course slightly change the underlying problem to some degree and make the formulation less clean. 

wM_std_VEC = [0.05]; %0.01 0.025 0.05 0.075 0.1
wPq_std_VEC = [3e-4]; % 6e-4];% 1.2e-3];% 2.4e-3];
wPqdot_std_VEC = [2.4e-3 ];%4.8e-3];% 9.6e-3];
addpath("~/casadi-3.6.5")
% addpath("./6muscle_control")
addpath("./forwardSim")
addpath("./Muscle_LMT_dM")
addpath("./MuscleModel")
addpath("./ArmModel")
addpath("./MusculoskeletalDynamics")
addpath("./Integrator")
addpath("./plotFunctions")
listing = dir();

wM_std = wM_std_VEC(1);
forceField = 0;

% saveName = ['result_time_0.8_' target '_forceField_' num2str(forceField) '_' num2str(wM_std) '_' num2str(wPq_std) '_' num2str(wPqdot_std) '.mat'];
result = optimization_6muscles(forceField,wM_std);

%%
plotResults(result);

%%
animate_trajectory(result);

%%
    subplot(2, 2, 1:4)
    for i = 1:numel(81)
        cla
        hold on 
        title("Planar Arm $t$ = " + num2str(i, '%.2f') + " s", 'Interpreter', 'latex')
        l1 = 1; l2 = 1.5;
        theta_shoulder_i = 0.01*i;
        theta_elbow_i = 0.01*i;
        link1 = [0 0; l1*cos(theta_shoulder_i) l1*sin(theta_shoulder_i)];
        link2 = [link1(2, :); link1(2, 1) + l2*cos(theta_shoulder_i + theta_elbow_i) link1(2, 2) + l2*sin(theta_shoulder_i + theta_elbow_i)];

        plot(link1(:, 1), link1(:, 2), 'b-', 'LineWidth', 4)
        plot(link2(:, 1), link2(:, 2), 'b-', 'LineWidth', 4)
        plot(link1(2, 1), link1(2, 2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
        plot(link2(2, 1), link2(2, 2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')

        xlim([-5, 5])
        ylim([-5, 5])
        hold off
        axis equal 
        axis([-5 5 -5 5])
        xlabel("$x$ position (m)", 'Interpreter', 'latex')
        ylabel("$y$ position (m)", 'Interpreter', 'latex')
        drawnow
    end
%%
