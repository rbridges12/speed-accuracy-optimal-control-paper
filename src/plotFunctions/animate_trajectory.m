function animate_trajectory(result)
    ts = result.time;
    q = result.q;
    P = result.Pmat;
    EEpos = result.EEPos;
    P_EEpos = result.P_EEPos;
    EEvel = result.EEVel;
    P_EEvel = result.P_EEVel;
    l1 = result.auxdata.l1;
    l2 = result.auxdata.l2;

    FPS = 200;
    t_anim = 0:1/FPS:ts(end);

    q_anim = [interp1(ts, q(:, 1), t_anim);
              interp1(ts, q(:, 2), t_anim)];
    EEpos_anim = [interp1(ts, EEpos(:, 1), t_anim);
                  interp1(ts, EEpos(:, 2), t_anim)];
    EEvel_anim = [interp1(ts, EEvel(:, 1), t_anim);
                  interp1(ts, EEvel(:, 2), t_anim)];
    sigma_1_anim = sqrt(interp1(ts, squeeze(P(1, 1, :)), t_anim));
    sigma_2_anim = sqrt(interp1(ts, squeeze(P(2, 2, :)), t_anim));
    pos_cov_anim = zeros(2, 2, numel(t_anim));
    pos_cov_anim(1, 1, :) = interp1(ts, P_EEpos(1, :), t_anim);
    pos_cov_anim(2, 2, :) = interp1(ts, P_EEpos(3, :), t_anim);
    cov_xy = interp1(ts, P_EEpos(2, :), t_anim);
    pos_cov_anim(1, 2, :) = cov_xy;
    pos_cov_anim(2, 1, :) = cov_xy;
    vel_cov_anim = zeros(2, 2, numel(t_anim));
    vel_cov_anim(1, 1, :) = interp1(ts, P_EEvel(1, :), t_anim);
    vel_cov_anim(2, 2, :) = interp1(ts, P_EEvel(3, :), t_anim);
    vel_cov_xy = interp1(ts, P_EEvel(2, :), t_anim);
    vel_cov_anim(1, 2, :) = vel_cov_xy;
    vel_cov_anim(2, 1, :) = vel_cov_xy;
    
    fig = figure;
    % video = VideoWriter('planar_arm_trajectory.mp4', 'MPEG-4')
    % video = VideoWriter('planar_arm_trajectory.avi')
    % video.FrameRate = FPS;
    % open(video)
    subplot(2, 2, 1:4)
    for i = 1:numel(t_anim)
        cla
        hold on 
        title("Arm Trajectory, final cost = " + num2str(result.final_cost, "%.2f") + ", t = " + num2str(t_anim(i), '%.2f') + " s")

        theta_shoulder_i = q_anim(1, i);
        theta_elbow_i = q_anim(2, i);

        % plot target
        draw_circle(result.EE_target, result.target_width, 'k--')

        plot(EEpos_anim(1, 1:i), EEpos_anim(2, 1:i), 'k-', 'LineWidth', 2)
        error_ellipse(pos_cov_anim(:, :, i), EEpos_anim(:, i), 0.95, "k:");

        EEx = EEpos_anim(1, i);
        EEy = EEpos_anim(2, i);
        quiver(EEx, EEy, -0.05*EEvel_anim(1, i), -0.05*EEvel_anim(2, i), 'r', 'LineWidth', 2)
        error_ellipse(vel_cov_anim(:, :, i).*0.05, [EEx; EEy], 0.95, "r:");

        sigma1_i = 2*sigma_1_anim(i);
        sigma2_i = 2*sigma_2_anim(i);
        draw_arm(theta_shoulder_i + sigma1_i, theta_elbow_i + sigma2_i, l1, l2, 2, 2, 'b--', 'co', 'ro', false)
        draw_arm(theta_shoulder_i - sigma1_i, theta_elbow_i - sigma2_i, l1, l2, 2, 2, 'b--', 'co', 'ro', false)

        plot([-0.05 0.05], [0 0], 'k-', 'LineWidth', 16)
        draw_arm(theta_shoulder_i, theta_elbow_i, l1, l2, 12, 8, 'b-', 'co', 'ro', true)

        legend("Target Area", "Mean EE Trajectory", "95% Confidence EE Position", "Mean EE Velocity", "95% Confidence EE Velocity", "95% Confidence Joint Angle Bounds")
        hold off
        axis equal 
        axis([-0.5 0.5 -0.1 0.8])
        xlabel("X Position (m)")
        ylabel("Y Position (m)")
        drawnow
        % frame = getframe(fig);
        % writeVideo(video, frame);
    end
    % close(video)
end

function draw_circle(center, radius, fmt)
    th = 0:pi/50:2*pi;
    xunit = radius * cos(th) + center(1);
    yunit = radius * sin(th) + center(2);
    plot(xunit, yunit, fmt)
end