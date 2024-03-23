function animate_trajectory(result)
    ts = result.time;
    q = result.q;
    EEpos = result.EEPos;
    P_EEpos = result.P_EEPos;
    l1 = result.auxdata.l1;
    l2 = result.auxdata.l2;

    FPS = 300;
    t_anim = 0:1/FPS:ts(end);

    q_anim = [interp1(ts, q(:, 1), t_anim);
              interp1(ts, q(:, 2), t_anim)];
    EEpos_anim = [interp1(ts, EEpos(:, 1), t_anim);
                  interp1(ts, EEpos(:, 2), t_anim)];
    cov_anim = zeros(2, 2, numel(t_anim));
    cov_anim(1, 1, :) = interp1(ts, P_EEpos(1, :), t_anim);
    cov_anim(2, 2, :) = interp1(ts, P_EEpos(3, :), t_anim);
    cov_xy = interp1(ts, P_EEpos(2, :), t_anim);
    cov_anim(1, 2, :) = cov_xy;
    cov_anim(2, 1, :) = cov_xy;
    
    fig = figure;
    % video = VideoWriter('planar_arm_trajectory.mp4', 'MPEG-4')
    % video = VideoWriter('planar_arm_trajectory.avi')
    % video.FrameRate = FPS;
    % open(video)
    subplot(2, 2, 1:4)
    for i = 1:numel(t_anim)
        cla
        hold on 
        title("Planar Arm $t$ = " + num2str(t_anim(i), '%.2f') + " s", 'Interpreter', 'latex')

        theta_shoulder_i = q_anim(1, i);
        theta_elbow_i = q_anim(2, i);
        link1 = [0 0; l1*cos(theta_shoulder_i) l1*sin(theta_shoulder_i)];
        link2 = [link1(2, :); link1(2, 1) + l2*cos(theta_shoulder_i + theta_elbow_i) link1(2, 2) + l2*sin(theta_shoulder_i + theta_elbow_i)];

        plot(EEpos_anim(1, 1:i), EEpos_anim(2, 1:i), 'k--', 'LineWidth', 2)
        error_ellipse(cov_anim(:, :, i), EEpos_anim(:, i), 0.95);

        plot([-0.05 0.05], [0 0], 'k-', 'LineWidth', 16)
        plot(link1(:, 1), link1(:, 2), 'b-', 'LineWidth', 12)
        plot(link2(:, 1), link2(:, 2), 'b-', 'LineWidth', 8)
        plot(0, 0, 'co', 'MarkerSize', 14, 'MarkerFaceColor', 'c')
        plot(link1(2, 1), link1(2, 2), 'co', 'MarkerSize', 12, 'MarkerFaceColor', 'c')
        plot(link2(2, 1), link2(2, 2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r')

        legend("End Effector Position", "95\% Confidence Ellipse", "", "", "", "", "", "", 'Interpreter', 'latex')
        hold off
        axis equal 
        axis([-0.5 0.5 -0.1 0.8])
        xlabel("$x$ position (m)", 'Interpreter', 'latex')
        ylabel("$y$ position (m)", 'Interpreter', 'latex')
        drawnow
        % frame = getframe(fig);
        % writeVideo(video, frame);
    end
    % close(video)
end