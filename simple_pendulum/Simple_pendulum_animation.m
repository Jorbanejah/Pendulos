%In the latest framework, the pendulum's equation is written
%as: x'' + g/L*sin(x) = 0. This equation can be solve using ode45 or
%other numerical methods such as implicit Euler method, Heun's method...
%Now, let's animated this simple concept: 

%Since x represents the angle formed with the horizontal line, we will
%refer to it as theta.
g = 9.81;
L=1;
t_init = [0 5]; %Time
initcondition = [pi/4; 0]; %Initial angle and angular velocity
initcondition1 = [pi/4; sqrt((2*g)/L)];
initcondition2 = [pi/4; sqrt(g/L)];
[t, theta] = ode45(@equation, t_init, initcondition); %Whithout inicial condition
[t1, theta1] = ode45(@equation, t_init, initcondition1);%Angular velocity
[t2, theta2] = ode45(@equation, t_init, initcondition2);

nFrames = length(t);
theta1_intercept = interp1(t1, theta1(:,1), linspace(0, 5, nFrames));
theta2_intercept = interp1(t2, theta2(:,1), linspace(0, 5, nFrames));

%Record simulation:

filename = 'Simple_pendulum.gif';
fig = figure(1);

for i = 1:length(t)
    clf; %Wipe the current figure clean before drawing the next frame
    
    x = L * sin(theta(i,1));
    y = - L * cos(theta(i,1));

    hold on
    
    %Graphics stile:

    set(gcf, 'Color', 'w'); % Background white
    color_bar = [0.2 0.2 0.2]; % Gray
    color_line = [0 0.4470 0.7410]; % Blue MATLAB   
    color_ball = [0.8500 0.3250 0.0980]; % Orange MATLAB
       
    sgtitle('Pendulum Comparison at Different Angular Velocities', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Subplot 1:
    subplot(3,1,1)
    plot([0 x], [0 y], '-', 'Color', color_bar, 'LineWidth', 2); hold on
    plot([-1.2 1.2], [0 0], '-', 'Color', color_line, 'LineWidth', 2);
    plot(x, y, 'o', 'MarkerSize', 10, 'MarkerFaceColor', color_ball, 'MarkerEdgeColor', 'k');
    plot(0, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k'); % Pivot
    title('No Angular Velocity', 'FontSize', 14, 'FontWeight', 'bold');
    axis equal; 
    axis([-1.2 1.2 -1.2 0.2]); 
    box off;

    % Subplot 2:
    subplot(3,1,2)
    theta1 = theta1_intercept(i);
    x1 = L * sin(theta1);
    y1 = -L * cos(theta1);
    plot([0 x1], [0 y1], '-', 'Color', color_bar, 'LineWidth', 2); hold on
    plot([-1.2 1.2], [0 0], '-', 'Color', color_line, 'LineWidth', 2);
    plot(x1, y1, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k');        
    plot(0, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
    title('Angular Velocity: (2g/L)', 'FontSize', 14, 'FontWeight', 'bold');
    axis equal; 
    axis([-1.2 1.2 -1.2 1.2]); 
    box off;

    % Subplot 3: Velocidad angular sqrt(g/L)
    subplot(3,1,3)
    theta2 = theta2_intercept(i);
    x2 = L * sin(theta2);
    y2 = -L * cos(theta2);
    plot([0 x2], [0 y2], '-', 'Color', color_bar, 'LineWidth', 2); hold on
    plot([-1.2 1.2], [0 0], '-', 'Color', color_line, 'LineWidth', 2);
    plot(x2, y2, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'w');
    plot(0, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
    title('Angular Velocity: √(g/L)', 'FontSize', 14, 'FontWeight', 'bold');

    axis equal; 
    axis([-1.2 1.2 -1.2 0.4]); 
    box off;
    
    drawnow;
    frame = getframe(fig);
    img = frame2im(frame);
    [A,map] = rgb2ind(img,256);

    if i == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.05);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.05);
    end
end
%First step is converted the second-order ODE into a first-order ODE
%through this function:
function dtheta = equation(t, theta)
    g = 9.81;
    L = 1;
    dtheta = zeros(2,1);
    dtheta(1) = theta(2);
    dtheta(2) = -(g/L)*sin(theta(1));
end