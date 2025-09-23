%In the latest framework, the pendulum's equation is written
%as: x'' + g/L*sin(x) = 0. This equation can be solve using ode45 or
%other numerical methods such as implicit Euler method, Heun's method...
%Now, let's animated this simple concept: 

%Since x represents the angle formed with the horizontal line, we will
%refer to it as theta.
g = 9.81;

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

L = 1;
for i = 1:length(t)
    clf; %Wipe the current figure clean before drawing the next frame
    grid on
    subplot(3,1,1);
    
    x = L * sin(theta(i,1));
    y = - L * cos(theta(i,1));

    hold on

    plot([0 x], [0 y], 'k-', 'LineWidth', 2); %Black line
    plot([-1.2 1.2], [0 0], 'b','LineWidth', 2);%Top blue line
    plot(x, y, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); %Ball
    title('Pendulum with no angular velocity')
    axis equal; axis([-1.2 1.2 -1.2 0.2]);

    subplot(3,1,2)

    x1 = L * sin(theta1_intercept(i));
    y1 = -L * cos(theta1_intercept(i));

    hold on

    plot([0 x1], [0 y1], 'k-', 'LineWidth', 2); %Black line
    plot([-1.2 1.2], [0 0], 'b','LineWidth', 2);%Top blue line
    plot(x1, y1, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'b'); %Ball
    title('Pendulum with sqrt(2g/L) angular velocity')
    axis equal; axis([-1.2 1.2 -1.2 1.2]);

    subplot(3,1,3);

    x2 = L * sin(theta2_intercept(i));
    y2 = -L * cos(theta2_intercept(i));

    hold on

    plot([0 x2], [0 y2], 'k-', 'LineWidth', 2); %Black line
    plot([-1.2 1.2], [0 0], 'b','LineWidth', 2);%Top blue line
    plot(x2, y2, 'ro', 'MarkerSize', 10, 'MarkerEdgeColor', 'g');%Ball
    title('Pendulum with sqrt(g/L) angular velocity')
    axis equal; axis([-1.2 1.2 -1.2 0.4]);

    drawnow;
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