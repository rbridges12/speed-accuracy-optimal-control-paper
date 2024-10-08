function plotNormalizedVelocity(result)
    EE_vel = result.EEVel;
    norm_vel = vecnorm(EE_vel,2,2);
    [max_vel, max_vel_i] = max(norm_vel);
    max_vel_time = result.time(max_vel_i)
    max_vel

    normalized_vel = norm_vel./max(norm_vel);
    normalized_time = result.time./max(result.time);

    figure
    plot(normalized_time, normalized_vel, 'LineWidth', 2)
    % plot(normalized_time, norm_vel, 'LineWidth', 2)
    title('Normalized Velocity of End Effector')
    xlabel('Normalized Time')
    ylabel('Normalized Velocity')
end