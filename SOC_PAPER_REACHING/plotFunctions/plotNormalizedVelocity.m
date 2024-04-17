function plotNormalizedVelocity(result)
    EE_vel = result.EEVel;
    norm_vel = vecnorm(EE_vel,2,2)

    figure
    plot(result.time, norm_vel, 'LineWidth', 2)
    title('Normalized Velocity of End Effector')
    xlabel('Time (s)')
    ylabel('Normalized Velocity (m/s)')
end