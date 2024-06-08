function draw_arm(theta_shoulder, theta_elbow, l1, l2, l1_size, l2_size, l_fmt, j_fmt, ee_fmt, print_angles)
    link1 = [0 0; l1*cos(theta_shoulder) l1*sin(theta_shoulder)];
    link2 = [link1(2, :); link1(2, 1) + l2*cos(theta_shoulder + theta_elbow) link1(2, 2) + l2*sin(theta_shoulder + theta_elbow)];

    plot(link1(:, 1), link1(:, 2), l_fmt, 'LineWidth', l1_size)
    plot(link2(:, 1), link2(:, 2), l_fmt, 'LineWidth', l2_size)
    plot(0, 0, j_fmt, 'MarkerSize', 14, 'MarkerFaceColor', 'c')
    plot(link1(2, 1), link1(2, 2), j_fmt, 'MarkerSize', l1_size, 'MarkerFaceColor', 'c')
    plot(link2(2, 1), link2(2, 2), ee_fmt, 'MarkerSize', l2_size, 'MarkerFaceColor', 'r')
    if print_angles
        rad2deg = 180/pi;
        x_offset = 0.02;
        y_offset = 0.005;
        text(link1(1, 1) - x_offset, link1(1, 2) - 0.04, "\theta_1 = " + num2str(rad2deg*theta_shoulder, '%.f') + "°", 'FontWeight', 'bold', 'FontSize', 12)
        text(link2(1, 1) + x_offset, link2(1, 2) + y_offset, "\theta_2 = " + num2str(rad2deg*theta_elbow, '%.f') + "°", 'FontWeight', 'bold', 'FontSize', 12)
    end
end