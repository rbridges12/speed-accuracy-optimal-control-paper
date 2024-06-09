function sigma_test(result)
    ts = result.time;
    q = result.q;
    X = result.X;
    u = result.e_ff;
    P = result.Pmat;
    EEpos = result.EEPos;
    P_EEpos = result.P_EEPos;
    EEvel = result.EEVel;
    P_EEvel = result.P_EEVel;
    l1 = result.auxdata.l1;
    l2 = result.auxdata.l2;
    nStates = result.auxdata.nStates;
    Q = result.Sigma_w;

    P_init = squeeze(P(:,:,1));
    x_init = X(:,1);
    P_init = blkdiag(P_init, Q);
    x_init = [x_init; zeros(size(Q,1),1)];
    kappa = 2;
    n = numel(x_init);
    dt = result.time(2);
    wM = zeros(result.auxdata.nMotorNoises,1);
    f = @(x, u, w) x + forwardMusculoskeletalDynamics_motorNoise(x, u, 0, w, result.auxdata) * dt;
    
    % sigma points around the reference point
    L = sqrt(n + kappa) * chol(P_init, 'lower');
    Y = x_init(:, ones(1, numel(x_init)));
    X = [x_init, Y + L, Y - L];
    w = zeros(2 * n + 1, 1);
    w(1) = kappa / (n + kappa);
    w(2:end) = 1 / (2*(n + kappa));

    Y = zeros(size(X));
    s_polar = (chol(P_init,'lower') * randn(n,1000) + x_init)';

    figure;
    for j = 1:numel(ts)
        mean = 0;
        for i = 1:2*n+1
            % Y(1:4,i) = f(X(1:4,i), u(j,:)', wM);
            Y(1:4,i) = f(X(1:4,i), u(j,:)', X(5:end,i));
            Y(5:end,i) = X(5:end,i);
            mean = mean + w(i) * Y(1:4,i);
        end
        % size(Y)
        Cov = (Y(1:4,:) - mean) * diag(w) * (Y(1:4,:) - mean)';
        % Cov_xy = (X - x) * diag(w) * (Y - mean)';

        subplot(2,1,1)
        cla
        hold on
        for i = 1:size(Y, 2)
            draw_arm(Y(1,i), Y(2,i), l1, l2, 2, 2, 'b-', 'co', 'ro', false)
        end
        axis equal 
        axis([-0.5 0.5 -0.1 0.8])
        xlabel("X Position (m)")
        ylabel("Y Position (m)")

        subplot(2,1,2)
        cla
        hold on
        fsize = 22; % font size
        set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
        set(groot, 'defaultLegendInterpreter','latex');
        Darkgrey = [.25 .25 .25];
        VermillionRed = [156,31,46]/255;
        DupontGray = [144,131,118]/255;

        % create confidence ellipse
        % first create points from a unit circle
        phi = (-pi:.01:pi)';
        circle = [cos(phi), sin(phi)];
        % Chi-squared 2-DOF 95% percent confidence (0.05): 5.991
        scale = sqrt(5.991);
        % apply the transformation and scale of the covariance
        % ellipse_polar = (scale * chol(P,'lower') * circle' + x)';
        mean_q = mean(1:2);
        Cov_q = Cov(1:2,1:2);
        ellipse_cartesian = (scale * chol(Cov_q,'lower') * circle' + mean_q)';

        % generate samples for both polar and cartesian coordinates
        s_cartesian = zeros(size(s_polar));
        for i = 1:size(s_polar,1)
            % s_cartesian(i,1:4) = f(s_polar(i,1:4)', u(j,:)', wM);
            s_cartesian(i,1:4) = f(s_polar(i,1:4)', u(j,:)', s_polar(i,5:end)');
            s_cartesian(i,5:end) = s_polar(i,5:end);
        end

        % plot in Cartesian coordinates
        % figure; hold on; grid on
        h = []; % plot handle
        h{1} = plot(s_cartesian(:,1), s_cartesian(:,2), '.', 'color', DupontGray);
        h{2} = plot(mean(1), mean(2), 'o', 'color', VermillionRed, 'markersize', 18);
        h{3} = plot(ellipse_cartesian(:,1), ellipse_cartesian(:,2), 'color', VermillionRed, 'linewidth', 3);
        h{4} = plot(Y(1,:), Y(2,:), '.', 'color', Darkgrey, 'markersize', 20);
        xlabel('$x=r\cos(\theta)$', 'Interpreter','latex'); 
        ylabel('$y=r\sin(\theta)$', 'Interpreter','latex');
        legend([h{1}, h{2}, h{3}, h{4}], 'Samples', 'Mean', '$95\%$ Confidence Ellipse', 'Sigma Points', 'location', 'north outside')
        % text(1.6, 1.8, '$\kappa = 2$', 'fontsize',fsize, 'Interpreter','latex')
        axis equal auto
        % set(gca,'fontsize',fsize)
        set(gca,'TickLabelInterpreter','latex')
        % axis([0.2 1.2 1.2 3])
        drawnow
        X = Y;
        s_polar = s_cartesian;
    end
end

%     for i = 1:numel(t_anim)
%         cla
%         hold on 
%         title("Arm Trajectory, final cost = " + num2str(result.final_cost, "%.2f") + ", t = " + num2str(t_anim(i), '%.2f') + " s")

%         theta_shoulder_i = q_anim(1, i);
%         theta_elbow_i = q_anim(2, i);

%         % plot target
%         draw_circle(result.EE_target, result.target_width, 'k--')

%         plot(EEpos_anim(1, 1:i), EEpos_anim(2, 1:i), 'k-', 'LineWidth', 2)
%         error_ellipse(pos_cov_anim(:, :, i), EEpos_anim(:, i), 0.95, "k:");

%         EEx = EEpos_anim(1, i);
%         EEy = EEpos_anim(2, i);
%         quiver(EEx, EEy, -0.05*EEvel_anim(1, i), -0.05*EEvel_anim(2, i), 'r', 'LineWidth', 2)
%         error_ellipse(vel_cov_anim(:, :, i).*0.05, [EEx; EEy], 0.95, "r:");

%         sigma1_i = 2*sigma_1_anim(i);
%         sigma2_i = 2*sigma_2_anim(i);
%         draw_arm(theta_shoulder_i + sigma1_i, theta_elbow_i + sigma2_i, l1, l2, 2, 2, 'b--', 'co', 'ro', false)
%         draw_arm(theta_shoulder_i - sigma1_i, theta_elbow_i - sigma2_i, l1, l2, 2, 2, 'b--', 'co', 'ro', false)

%         plot([-0.05 0.05], [0 0], 'k-', 'LineWidth', 16)
%         draw_arm(theta_shoulder_i, theta_elbow_i, l1, l2, 12, 8, 'b-', 'co', 'ro', true)

%         legend("Target Area", "Mean EE Trajectory", "95% Confidence EE Position", "Mean EE Velocity", "95% Confidence EE Velocity", "95% Confidence Joint Angle Bounds")
%         hold off
%         axis equal 
%         axis([-0.5 0.5 -0.1 0.8])
%         xlabel("X Position (m)")
%         ylabel("Y Position (m)")
%         drawnow
%         % frame = getframe(fig);
%         % writeVideo(video, frame);
%     end
%     % close(video)
% end


%     %% 2D visualization
%     fsize = 22; % font size
%     set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
%     set(groot, 'defaultLegendInterpreter','latex');
%     green = [0.2980 .6 0];
%     % crimson = [220,20,60]/255; 
%     % darkblue = [0 .2 .4];
%     Darkgrey = [.25 .25 .25];
%     % darkgrey = [.35 .35 .35];
%     % lightgrey = [.7 .7 .7];
%     % Lightgrey = [.9 .9 .9];
%     VermillionRed = [156,31,46]/255;
%     DupontGray = [144,131,118]/255;
%     % Azure = [53, 112, 188]/255;
%     % purple = [178, 102, 255]/255;
%     % orange = [255,110,0]/255;

%     % create confidence ellipse
%     % first create points from a unit circle
%     phi = (-pi:.01:pi)';
%     circle = [cos(phi), sin(phi)];
%     % Chi-squared 2-DOF 95% percent confidence (0.05): 5.991
%     scale = sqrt(5.991);
%     % apply the transformation and scale of the covariance
%     ellipse_polar = (scale * chol(P,'lower') * circle' + x)';
%     ellipse_cartesian = (scale * chol(Cov,'lower') * circle' + mean)';

%     % generate samples for both polar and cartesian coordinates
%     s_polar = (chol(P,'lower') * randn(2,1000) + x)';
%     s_cartesian = zeros(size(s_polar));
%     for i = 1:size(s_polar,1)
%         s_cartesian(i,:) = f(s_polar(i,:));
%     end

%     % plot in polar coordinates
%     figure; hold on; grid on
%     h = []; % plot handle
%     h{1} = plot(s_polar(:,1), s_polar(:,2), '.', 'color', DupontGray);
%     h{2} = plot(x(1), x(2), 'o', 'color', VermillionRed, 'markersize', 18);
%     h{3} = plot(ellipse_polar(:,1), ellipse_polar(:,2), 'color', VermillionRed, 'linewidth', 3);
%     h{4} = plot(X(1,:), X(2,:), '.', 'color', Darkgrey, 'markersize', 32);
%     xlabel('$r$', 'Interpreter','latex'); 
%     ylabel('$\theta$', 'Interpreter','latex');
%     legend([h{1}, h{2}, h{3}, h{4}], 'Samples', 'Mean', '$95\%$ Confidence Ellipse', 'Sigma Points', 'location', 'north outside')
%     text(1.75, 1.6, '$\kappa = 2$', 'fontsize',fsize, 'Interpreter','latex')
%     axis equal auto
%     set(gca,'fontsize',fsize)
%     set(gca,'TickLabelInterpreter','latex')
%     % figuresize(21,21,'cm')
%     % print -opengl -dpng -r600 ut_example_polar.png

%     % plot in Cartesian coordinates
%     figure; hold on; grid on
%     h = []; % plot handle
%     h{1} = plot(s_cartesian(:,1), s_cartesian(:,2), '.', 'color', DupontGray);
%     h{2} = plot(mean(1), mean(2), 'o', 'color', VermillionRed, 'markersize', 18);
%     h{3} = plot(ellipse_cartesian(:,1), ellipse_cartesian(:,2), 'color', VermillionRed, 'linewidth', 3);
%     h{4} = plot(Y(1,:), Y(2,:), '.', 'color', Darkgrey, 'markersize', 32);
%     xlabel('$x=r\cos(\theta)$', 'Interpreter','latex'); 
%     ylabel('$y=r\sin(\theta)$', 'Interpreter','latex');
%     legend([h{1}, h{2}, h{3}, h{4}], 'Samples', 'Mean', '$95\%$ Confidence Ellipse', 'Sigma Points', 'location', 'north outside')
%     text(1.6, 1.8, '$\kappa = 2$', 'fontsize',fsize, 'Interpreter','latex')
%     axis equal auto
%     set(gca,'fontsize',fsize)
%     set(gca,'TickLabelInterpreter','latex')
%     % figuresize(21,21,'cm')
% end