function sigma_test(result)
    ts = result.time;
    q = result.q;
    X = result.X;
    u = result.e_ff;
    P = result.Pmat;
    EEpos = result.EEPos';
    EEpos_vars = result.P_EEPos;
    P_EEpos = zeros(2, 2, numel(ts));
    P_EEpos(1, 1, :) = EEpos_vars(1, :);
    P_EEpos(2, 2, :) = EEpos_vars(3, :);
    P_EEpos(1, 2, :) = EEpos_vars(2, :);
    P_EEpos(2, 1, :) = EEpos_vars(2, :);
    l1 = result.auxdata.l1;
    l2 = result.auxdata.l2;
    nStates = result.auxdata.nStates;
    Q = result.Sigma_w;

    P_init = squeeze(P(:,:,1));
    x_init = X(:,1);
    P_init = blkdiag(P_init, Q);
    x_init = [x_init; zeros(size(Q,1),1)];
    x_sigma = x_init;
    P_sigma = P_init;
    kappa = 2;
    n_particles = 1000;
    n = numel(x_init);
    dt = result.time(2);
    wM = zeros(result.auxdata.nMotorNoises,1);
    f = @(x, u, w) x + forwardMusculoskeletalDynamics_motorNoise(x, u, 0, w, result.auxdata) * dt;
    
    particles = (chol(P_init,'lower') * randn(n,n_particles) + x_init);
    s_EEpos = zeros(2, n_particles);

    figure
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'defaultLegendInterpreter','latex');
    alpha = 0.7;
    Darkgrey = [.25 .25 .25];
    DarkBlue = [0 .2 .4];
    DarkBlue = [DarkBlue, alpha];
    Azure = [53, 112, 188]/255;
    Azure = [Azure, alpha];
    VermillionRed = [156,31,46]/255;
    VermillionRed = [VermillionRed, alpha];
    DupontGray = [144,131,118]/255;

    for j = 1:numel(ts)
        % extract linearized data
        X_j = X(:,j);
        P_j = squeeze(P(:,:,j));
        EEpos_j = EEpos(:,j);
        P_EEpos_j = P_EEpos(:,:,j);

        % propagate sigma points
        L = sqrt(n + kappa) * chol(P_sigma, 'lower');
        Y = x_sigma(:, ones(1, numel(x_sigma)));
        sigma_points = [x_sigma, Y + L, Y - L];
        w = zeros(2 * n + 1, 1);
        w(1) = kappa / (n + kappa);
        w(2:end) = 1 / (2*(n + kappa));
        sigma_points_new = zeros(size(sigma_points));
        u_EEpos = zeros(2, numel(w));
        sigma_mean = 0;
        EEpos_mean = 0;
        for i = 1:2*n+1
            % sigma_points_new(1:4,i) = f(sigma_points(1:4,i), u(j,:)', wM);
            sigma_points_new(1:4,i) = f(sigma_points(1:4,i), u(j,:)', sigma_points(5:end,i));
            sigma_points_new(5:end,i) = sigma_points(5:end,i);
            sigma_mean = sigma_mean + w(i) * sigma_points_new(1:4,i);

            u_EEpos(:,i) = EndEffectorPos(sigma_points_new(1:2,i),result.auxdata);
            EEpos_mean = EEpos_mean + w(i) * u_EEpos(:,i);
        end
        P_sigma(1:4, 1:4) = (sigma_points_new(1:4,:) - sigma_mean) * diag(w) * (sigma_points_new(1:4,:) - sigma_mean)';
        x_sigma(1:4) = sigma_mean;
        EEpos_cov = (u_EEpos - EEpos_mean) * diag(w) * (u_EEpos - EEpos_mean)';

        % propagate Monte Carlo simulation
        particles_new = zeros(size(particles));
        for i = 1:n_particles
            particles_new(1:4,i) = f(particles(1:4,i), u(j,:)', particles(5:end,i));
            particles_new(5:end,i) = mvnrnd(wM, Q, 1);
            s_EEpos(:,i) = EndEffectorPos(particles_new(1:2,i),result.auxdata);
        end
        s_mean = mean(particles_new(1:4,:),2);
        s_cov = cov(particles_new(1:4,:)');
        s_EEpos_mean = mean(s_EEpos,2);
        s_EEpos_cov = cov(s_EEpos');

        % plot arm
        t = tiledlayout(2, 2);
        nexttile([1 2]); cla; hold on
        title("Arm Trajectory")

        % plot target
        draw_circle(result.EE_target, result.target_width, 'k--')

        % draw sigma point arms
        % for i = 1:size(sigma_points_new, 2)
        %     weight = w(i) * 15;
        %     draw_arm(sigma_points_new(1,i), sigma_points_new(2,i), l1, l2, weight, weight, 'b--', 'co', 'ro', false)
        % end

        % draw mean arms
        % draw_arm(X_j(1), X_j(2), l1, l2, 1, 1, 'b-', 'co', 'ro', true)
        % draw_arm(sigma_mean(1), sigma_mean(2), l1, l2, 1, 1, 'r-', 'co', 'ro', true)
        plot([-0.05 0.05], [0 0], 'k-', 'LineWidth', 16)
        draw_arm(s_mean(1), s_mean(2), l1, l2, 12, 8, 'b-', 'co', 'ro', true)

        % formatting
        legend("Target")
        axis equal 
        axis([-0.5 0.5 -0.1 0.8])
        xlabel("X Position (m)")
        ylabel("Y Position (m)")

        % plot Monte Carlo simulation
        % nexttile; cla; hold on
        % title("Monte Carlo Simulation")

        % % plot target
        % draw_circle(result.EE_target, result.target_width, 'k--')


        % for i = 1:25:n_particles
        %     weight = 0.8;
        %     draw_arm(particles_new(1,i), particles_new(2,i), l1, l2, weight, weight, 'b--', 'co', 'ro', false)
        % end

        % % draw mean arm
        % draw_arm(s_mean(1), s_mean(2), l1, l2, 1, 1, 'r-', 'co', 'ro', true)

        % % formatting
        % axis equal 
        % axis([-0.5 0.5 -0.1 0.8])
        % xlabel("X Position (m)")
        % ylabel("Y Position (m)")

        % plot state space distributions
        nexttile; cla; hold on
        title("State Space Distribution")

        ellipse_linearized = confidence_ellipse(X_j(1:2), P_j(1:2,1:2));
        ellipse_unscented = confidence_ellipse(sigma_mean(1:2), P_sigma(1:2,1:2));
        ellipse_monte_carlo = confidence_ellipse(s_mean(1:2), s_cov(1:2,1:2));

        plot(particles_new(1,:), particles_new(2,:), '.', 'color', DupontGray);
        plot(sigma_points_new(1,:), sigma_points_new(2,:), '.', 'color', Darkgrey, 'markersize', 20);

        % linearized
        plot(X_j(1), X_j(2), 'o', 'color', Azure, 'markersize', 18, 'LineWidth', 1.5);
        plot(ellipse_linearized(:,1), ellipse_linearized(:,2), 'color', Azure, 'linewidth', 3);
        
        % unscented
        plot(sigma_mean(1), sigma_mean(2), 'o', 'color', VermillionRed, 'markersize', 18, 'LineWidth', 1.5);
        plot(ellipse_unscented(:,1), ellipse_unscented(:,2), 'color', VermillionRed, 'linewidth', 3);

        % monte carlo
        plot(s_mean(1), s_mean(2), 'o', 'color', DarkBlue, 'markersize', 18, 'LineWidth', 1.5);
        plot(ellipse_monte_carlo(:,1), ellipse_monte_carlo(:,2), 'color', DarkBlue, 'linewidth', 3);

        % formatting
        ylabel('$\theta_2$', 'Interpreter','latex');
        xlabel('$\theta_1$', 'Interpreter','latex'); 
        legend("Particles", "Sigma Points", "Linearized Mean", "Linearized Covariance", "Unscented Mean", "Unscented Covariance", "Monte Carlo Mean", "Monte Carlo Covariance", 'location', 'east outside')
        axis equal auto
        set(gca,'TickLabelInterpreter','latex')

        % plot end effector space distributions
        nexttile; cla; hold on
        title("End Effector Distribution")

        ellipse_EE_linearized = confidence_ellipse(EEpos_j, P_EEpos_j);
        ellipse_EE_unscented = confidence_ellipse(EEpos_mean, EEpos_cov);
        ellipse_EE_monte_carlo = confidence_ellipse(s_EEpos_mean, s_EEpos_cov);

        plot(s_EEpos(1,:), s_EEpos(2,:), '.', 'color', DupontGray);
        plot(u_EEpos(1,:), u_EEpos(2,:), '.', 'color', Darkgrey, 'markersize', 20);

        plot(EEpos_j(1), EEpos_j(2), 'o', 'color', Azure, 'markersize', 18, 'LineWidth', 1.5);
        plot(ellipse_EE_linearized(:,1), ellipse_EE_linearized(:,2), 'color', Azure, 'linewidth', 3);

        plot(EEpos_mean(1), EEpos_mean(2), 'o', 'color', VermillionRed, 'markersize', 18, 'LineWidth', 1.5);
        plot(ellipse_EE_unscented(:,1), ellipse_EE_unscented(:,2), 'color', VermillionRed, 'linewidth', 3);

        plot(s_EEpos_mean(1), s_EEpos_mean(2), 'o', 'color', DarkBlue, 'markersize', 18, 'LineWidth', 1.5)
        plot(ellipse_EE_monte_carlo(:,1), ellipse_EE_monte_carlo(:,2), 'color', DarkBlue, 'linewidth', 3);

        axis equal auto

        drawnow
        % sigma_points = sigma_points_new;
        particles = particles_new;
    end
end

function draw_circle(center, radius, fmt)
    th = 0:pi/50:2*pi;
    xunit = radius * cos(th) + center(1);
    yunit = radius * sin(th) + center(2);
    plot(xunit, yunit, fmt)
end

function ellipse = confidence_ellipse(mean, cov)
    phi = (-pi:.01:pi)';
    circle = [cos(phi), sin(phi)];
    % Chi-squared 2-DOF 95% percent confidence (0.05): 5.991
    scale = sqrt(5.991);
    ellipse = (scale * chol(cov,'lower') * circle' + mean)';
end