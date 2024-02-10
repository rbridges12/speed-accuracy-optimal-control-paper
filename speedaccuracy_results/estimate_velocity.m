%pos_data = Hs{9};

function vel_data = estimate_velocity(pos_data, time)

vel_data = zeros(size(pos_data,1),1);

for i = 2 : size(pos_data,1)
    dt = time(i) - time(i-1);
    vel_data(i) = (norm(pos_data(i,1:3) - pos_data(i-1,1:3)))/dt;  
    %vel_data(i) = ((norm(pos_data(i,1) - pos_data(i-1,1))) ...
    %    + (norm(pos_data(i,3) - pos_data(i-1,3))))/dt; 
end

vel_data = smooth(vel_data);

%{
width = 3;     % Width in inches
height = 3;    % Height in inches
alw = 2.75;    % AxesLineWidth
fsz = 24;      % Fontsize
lw = 4.5;      % LineWidthopen
msz = 8;       % MarkerSize

plot(time(1:end-1), vel_data, 'LineWidth', lw);
set(gca, 'FontSize', fsz, 'LineWidth', lw)
xlabel('Time')
ylabel('Velocity')

%saveas(gcf, 'vel4', 'fig');
%saveas(gcf, 'vel4', 'epsc');
%}