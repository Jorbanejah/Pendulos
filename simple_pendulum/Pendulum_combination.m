clc;clear
l = 1; M = 2; m = 1; g = 9.81; k = 2;
b = 0.2; %Cart dimension
t_init = [0 10];
conditional = [0; 0; pi/4; 0]; %Initial position, inicial velocity, angular position, angular velocity.
[t, z] = ode45(@(t, y) pendulum(t, y, g, l, m, M, k), t_init, conditional);

v = VideoWriter('Combination_simulation.gif');
v.FrameRate = 20; 
open(v);

nFrames = length(t);
for i = 1:nFrames
    clf;
    hold on 
    cart_x = z(i,1); %Cart
    theta = z(i,3); %Angle formed by the pendulum

    %We write the coordenates of the system:
    
    pend_x = cart_x + l * sin(theta);
    pend_y = 0.2 - l * cos(theta);  % asuming pivot height is 0.2
    
    %Cart: define the four points and then draw inside with fill

    cx = [cart_x-b, cart_x+b, cart_x+b, cart_x-b];
    cy = [0,        0,        2*b,      2*b];
    fill(cx, cy, [0.2 0.6 1], 'EdgeColor','k','LineWidth',2);
    
    %Pendulum

    plot([cart_x pend_x], [b pend_y], 'k-', 'LineWidth', 2);
    plot(pend_x, pend_y, 'ro', 'MarkerSize', 10,  'MarkerFaceColor', 'r');
    plot([-1.2 1.2], [0 0], 'k', 'LineWidth', 2);

   % Spring effect: use sine function to generate a coil-shaped spring
    nCoils = 6;
    L0 = 2;
    spring_y = b;
    Lcurr  = abs(cart_x) + 0.05;      % Current length
    tau    = linspace(0,1,200);
    sx     = tau * cart_x;           
    sy     = spring_y + 0.02*sin(2 * pi * nCoils * (tau * L0 / Lcurr));
    plot(sx, sy, 'g', 'LineWidth', 2);

    axis equal
    xlim([-1.2 1.2]);
    ylim([-1.2 1.2]);

    hold off
    drawnow;
    frame = getframe(gcf);
    writeVideo(v, frame);
end
close(v)
function dy = pendulum(t, y, g, L, m, M, k)

    x = y(1);
    dxdt = y(2);
    theta = y(3);
    dtheta = y(4);

    A = [(m + M) (M * L * cos(theta));
        (m * L * cos(theta)) (M * L^2)];

    B = [(M * L * dtheta^2 * sin(theta) - k * x);
        (x * M * L * dtheta * sin(theta) - M * L * x * dtheta^2 *sin(theta) - M * g * L * sin(theta))];

    sol = A \ B;
    dxdt_sol = sol(1);
    dtheta_sol = sol(2);

    dy = zeros(4,1);
    dy(1) = dxdt ;
    dy(2) = dxdt_sol;
    dy(3) = dtheta;
    dy(4) = dtheta_sol;
end